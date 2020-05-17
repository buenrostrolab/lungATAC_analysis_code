### Author: Vinay Kartha
### Contact: <vinay_kartha@g.harvard.edu>
### Affiliation: Buenrostro Lab, Department of Stem Cell and Regerative Biology, Harvard University

library(dplyr)
library(doParallel)
library(foreach)
library(scales)
library(genefilter)
library(ggplot2)

   
# Function to get TF short names frrom chromVARmotifs package attribute names
extractTFnames <-  function(motifIDs){
    sapply(strsplit(sapply(strsplit(motifIDs,"_LINE.",fixed=FALSE),"[[",2),"_",fixed=FALSE),"[[",2)
  }

# For a given motif Z score matrix, return a sparse binary matrix 
# Indicating high (1) vs low (0) median splits of cells based on Z score
.binarizeZMat <- function(Z){
  (Z > matrixStats::rowMedians(Z))*1
}


sigVolcano <- function(motifList,
                       motif=NULL,
                       FDR=0.000001){
  if(is.null(motif))
    stop("Must specify a valid motif ID to generate volcano plot for ..\n")
  
  d <- motifList[[motif]]
  
  # FDR from p-val
  d$FDR <- p.adjust(d$p.value,method="fdr")
  if(any(is.na(d$FDR))){
    cat(sum(is.na(d$FDR))," NA FDR values from t-test detected ..\n",sep = " ")
    cat("Removing these peaks prior to visualization..\n")
    d <- d[!is.na(d$FDR),]
  }
  d$Sig <- factor(ifelse(d$FDR < FDR,"Yes","No"),levels=c("Yes","No"))
  numSig <- sum(d$FDR < FDR)
  d$logFDR <- -log10(d$FDR)
  
  num0 <- log10(FDR) # Just to display the -ve 10th power integer instead of scientific format
  ggplot(d,aes(x=dm,y=logFDR,color=Sig)) + 
    ggrastr::geom_point_rast(size=0.2,alpha=0.8) + 
    theme_bw() + scale_color_manual(values=c("mediumpurple4","lightsteelblue")) + 
    labs(title=BuenRTools::extractTFNames(motif),
         subtitle=substitute(paste("Number of peaks with FDR < ",10^num0,": ",numSig,sep = ""),list(num0=num0,numSig=numSig)),
         x="Difference in mean accesibility (high - low)",y="-log10(FDR)") + 
    theme(plot.title=element_text(face="bold.italic",hjust=0.5),plot.subtitle=element_text(hjust=0.5)) + 
    guides(colour = guide_legend(override.aes = list(size=3),ncol=1)) + scale_y_continuous(limits = c(0,max(d$logFDR+3)),expand = c(0,0))
}

# NOTE: USE THIS TEST WHEN BINNING CELLS AS HIGH/LOW AND DOING T-TEST OF ASSOCIATION COMPARING COUNTS BETWEEN GROUPS
ttestPeaksMatPar <- function(Zscores, # Matrix of motifs x cells (Z scores)
                             scSE, # Matrix of peak x cells counts (scATAC-seq reads)
                             binarizeMat=FALSE, # Whether or not to binarize matrix?
                             normalizeMat=TRUE, # Whether or not to normalize peak counts per cell by mean
                             chunkSize=20000, # Number of peaks to test at once (in series or in parallel)
                             ncores=1, # Cores to use if parallelizing
                             byMotifs=TRUE){ # Return list of motifs, with all peak tests per motif
  stopifnot(all.equal(colnames(Zscores),colnames(scSE)))
  
  if(normalizeMat){
    cat("Normalizing matrix of counts based on mean reads in peaks..\n")
    scSE <- centerCounts(scSE)
  } else if(binarizeMat){
    cat("Binarizing matrix of counts ..\n")
    assay(scSE) <- (assay(scSE) > 0)*1
  } 
  
  cat("Getting high vs low groupings for cells based on median Z scores per motif..\n")
  Zgroups <- .binarizeZMat(Zscores)
  
  cat("Testing peaks sequentially within chunks of size ",
      chunkSize, " ..\n\n")
  
  if(ncores > 1)
    cat("Running chunks in parallel using ",ncores, "cores ..\n")
  
  starts <- seq(1,nrow(scSE),chunkSize)
  ends <- starts + chunkSize -1
  ends[length(ends)] <- nrow(scSE)
  
  chunkList <- mapply(c, starts, ends, SIMPLIFY = FALSE)
  
  time_elapsed <- Sys.time()
  
  cat("Running peakwise motif-peak accessibility association testing ..\n")
  peakList <- parallel::mclapply(chunkList,FUN = function(chunk,SE=scSE,Z=Zgroups){
    
    SE.chunk <- SE[chunk[1]:chunk[2],] # Subset chunk of peaks (rows)
    
    # Run row t-test on chunk with all motifs
    lapply(1:nrow(Z),FUN = function(motif,counts=as.matrix(assay(SE.chunk)),Zstatus=Z){
      m <- rownames(Zstatus)[motif]
      d <- genefilter::rowttests(counts,fac = factor(Zstatus[motif,])) # This returns dm (difference in mean) 0 (low group) - 1 (high group)
      d$motif <- m
      d
    })
    
  },mc.cores = ncores)
  
  
  cat("Merging results ..\n")
  merged.tab <- dplyr::bind_rows(unlist(peakList,recursive=FALSE))
  cat("Finished!\n\n")
  
  time_elapsed <- Sys.time() - time_elapsed
  cat(paste("\nTime Elapsed: ",time_elapsed, units(time_elapsed), "\n"))
  
  
  
  if(byMotifs)
    return(split.data.frame(merged.tab,
                            f=factor(merged.tab$motif,levels=unique(merged.tab$motif))))
  
  return(merged.tab)
}
