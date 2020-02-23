library(irlba)
library(foreach)
library(doParallel)

EmpiricalP <- function(obsval, nullval) {
  return(sum(nullval > obsval) / length(nullval))
}

MatrixRowShuffle <- function(x) {
  x2 <- x
  x2 <- t(x = x)
  ind <- order(c(col(x = x2)), runif(n = length(x = x2)))
  x2 <- matrix(
    data = x2[ind],
    nrow = nrow(x = x),
    ncol = ncol(x = x),
    byrow = TRUE
  )
  return(x2)
}


# Main PCA function
runPCA <- function(mat,
                   nPCs=20,
                   scale=FALSE,
                   center=FALSE){

  pc <- irlba::irlba((scale(mat,center = center,scale = scale)),nv = nPCs)

  # Matrix of gene/motif loadings (weighted by variance)
  motif.loadings <- pc$u %*% diag(pc$d)
  rownames(motif.loadings) <- rownames(mat) # This assumes the rows are motif IDs
  colnames(motif.loadings) <- paste0("PC",1:nPCs)
  # Matrix of PCA embeddings
  cell.embeddings <- pc$v

  return(list(motifLoadings=motif.loadings,
              cellEmbeddings=cell.embeddings))

}

# Jackstraw version of regular PCA run

JackPCA  <- function(mat,
                     propFeatures=0.25, # Proportion of total feature space to use when running jackstraw PCA
                     useUptoPC=20,
                     seed=123,
                     scale=FALSE,
                     center=FALSE){
  set.seed(seed)
  rand.features <- sample(
    x = rownames(mat),
    size = nrow(mat) * propFeatures
  )

  # make sure that rand.genes is at least 3
  if (length(x = rand.features) < 3){
    warning("Too few features selected during jackstraw subsampling procedure ..\n Consider testing more features, or increasing the propFeatures percentage.")
    rand.features <- sample(x = rownames(x = mat), size = 3)
  }


  # Randomly permute data in matrix (independent row shuffling) to disrupt data
  mat.rand <- MatrixRowShuffle(mat[rand.features,])
  rownames(mat.rand) <- rand.features

  pca.jack <- runPCA(mat.rand,nPCs = useUptoPC,scale = scale,center = center)

  jack.loadings <- pca.jack$motifLoadings[rand.features,1:useUptoPC]
  return(jack.loadings)
}


# Wrapper for running jackstraw version of PCA for a number of iterations. Main workflow and functions adapted from the Seurat R package.
jackstrawMotifs <- function(mat,
                            motifLoadings=NULL,
                            nPCs=20,
                            propUse=0.25,
                            nIterations=100,
                            do.par=FALSE,
                            num.cores=1,
                            display.progress=TRUE,
                            scale=FALSE,
                            center=FALSE){

  if(is.null(motifLoadings)){
    cat("Running first pass PCA on data ..\n")
    motifLoadings <- runPCA(mat = mat,nPCs = nPCs,scale = scale,center = center)$motifLoadings
  }

  # input checking for parallel options
  if (do.par) {
    if (num.cores == 1) {
      num.cores <- parallel::detectCores() / 2
      warning(paste0("do.par set to TRUE but num.cores set to 1. Setting num.cores to ", num.cores, "."))
    } else if (num.cores > detectCores()) {
      num.cores <- detectCores() - 1
      warning(paste0("num.cores set greater than number of available cores(", detectCores(), "). Setting num.cores to ", num.cores, "."))
    }
  } else if (num.cores != 1) {
    num.cores <- 1
    warning("For parallel processing, please set do.par to TRUE.")
  }


  cl <- parallel::makeCluster(num.cores)

  doSNOW::registerDoSNOW(cl)

  cat("Running jackstraw PCA to determine significant PCs ..\n")
  if (display.progress) {
    time_elapsed <- Sys.time()
  }

  opts <- list()
  if (display.progress) {
    # define progress bar function
    pb <- txtProgressBar(min = 0, max = nIterations, style = 3)
    progress <- function(n) setTxtProgressBar(pb, n)
    opts <- list(progress = progress)
    time_elapsed <- Sys.time()
  }

  jack.pcVals.raw <- foreach(
    x = 1:nIterations,
    .options.snow = opts,
    .export = c('JackPCA','EmpiricalP','MatrixRowShuffle','runPCA')
  ) %dopar% {
    JackPCA(mat=mat,
      propFeatures = propUse,
      useUptoPC = nPCs,
      seed=x,
      scale = scale,
      center = center
    )
  }

  if (display.progress) {
    time_elapsed <- Sys.time() - time_elapsed
    cat(paste("\nTime Elapsed: ",time_elapsed, units(time_elapsed), "\n"))
    close(pb)
  }

  parallel::stopCluster(cl)

  jack.pcVals <- sapply( # This will be a final matrix of concatenated scores from each iteration (rows) by PCs
    X = 1:nPCs,
    FUN = function(x) {
      return(as.numeric(unlist(lapply(
        X = 1:nIterations,
        FUN = function(y) {
          return(jack.pcVals.raw[[y]][, x]) # Extract motif loadings for each PC per iteration
        }
      ))))
    }
  )

  jackStraw.PCmat <- as.matrix(jack.pcVals)
  jackStraw.empP <- as.matrix(
    sapply(
      X = 1:nPCs,
      FUN = function(x) {
        return(unlist(x = lapply(
          X = abs(motifLoadings[, x]),
          FUN = EmpiricalP,
          nullval = abs(jack.pcVals[,x])
        )))
      }
    )
  )
  colnames(jackStraw.empP) <- paste0("PC", 1:ncol(jackStraw.empP))

  return(list("jackStrawPCScores"=jackStraw.PCmat,
       "jackStrawPCEmpPvals"=jackStraw.empP))

}


# This function takes the min emp pvalue estimated from jackstraw PCs (across all specified PCs per motif) and returns motifs passing a pval cutoff
getSigMotifs <- function(JSpvals, # Jackstraw pvals from jackstraw PCA on data
                         pcs.use=1:10,
                         pval.cutoff=0.1,
                         max.per.pc=NULL
                         ) {

  if (length(pcs.use) == 1) {
    pvals.min <- JSpvals[, pcs.use] # Usually just the first PC
  }
  if (length(pcs.use) > 1) {
    pvals.min <- apply(X = JSpvals[, pcs.use], MARGIN = 1, FUN = min)
  }

  names(pvals.min) <- rownames(JSpvals)
  sig.motifs <- names(pvals.min)[pvals.min < pval.cutoff]
  return(sig.motifs)
}

runTest <- FALSE
if(runTest){
chromdev <- readRDS("./motifDev_bagLeaders.rds")

Zscores <- scale(deviationScores(chromdev)[-1,],center = TRUE,scale = TRUE)

# Run jackstraw PCA for scaled/centered Z-scores of bag leaders
JSM <- jackstrawMotifs(mat = Zscores,
                       propUse = 0.2,
                       nPCs = 20,
                       nIterations = 1000)

sigmotifs <- getSigMotifs(JSpvals = JSM$jackStrawPCEmpPvals,pcs.use = 1:10,pval.cutoff = 0.05)
}