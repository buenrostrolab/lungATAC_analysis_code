# Source module functions
source("<path_to_git>/R/moduleDiffPeaks.R")
# Source jackstraw functions
source("<path_to_git>/R/jackStrawPCA.R.R")

setwd("<data_analysis_folder>")

# Load SE object
SE <- readRDS("./atac.se.rds")

# Load motif deviations object
devMotif <- readRDS("./devMotif.rds")

# Filter cells used in module analysis (remove normal cells)
cellsToKeep <- readRDS("./lungMetCells.rds")

SE <- SE[,cellsToKeep]
devMotif <- devMotif[,cellsToKeep]

# Bag motifs based on sequence similarity
devBagged <- bagDeviations(devMotif,cor = 0.8,organism = "human")

# Fetch scaled motif Z-scores
Zscaled <- scale(deviationScores(devBagged),center = TRUE,scale = TRUE)
stopifnot(sum(is.na(Zscaled)==0))

# Run jackstraw to keep sig variable motifs
JSM <- jackstrawMotifs(mat = Zscaled,
                       propUse = 0.20,
                       nPCs = 20,
                       nIterations = 1000,num.cores = 4,do.par = TRUE)

# Fetch sig variable motifs
jackSigmotifs <- getSigMotifs(JSpvals = JSM$jackStrawPCEmpPvals,pcs.use = 1:10,pval.cutoff = 0.05)

# Only use jackstraw PC sig motifs, and filtered cells
Z <- deviationScores(devBagged)[jackSigmotifs,]

saveRDS(Z,"./chromVAR/jackSigMelanomaMotifDev.rds")

# Do module analysis on bagged filtered motifs/cells
cat("Starting module analyses ..\n")
peakTFMotifList <- ttestPeaksMatParNew(Zscores = Z,
                                      scSE = SE,
                                      binarizeMat = FALSE,
                                      normalizeMat = TRUE,
                                      ncores = 6,
                                      byMotifs = TRUE)
cat("Finished!\n")

# REVERSE dm sign (making it high - low, instead of default low - high)
peakTFMotifList <- lapply(peakTFMotifList,function(m) { m$dm <- -(m$dm); return(m)})

# # Add FDR per motif test to adjust p-val
FDR <- 1e-06
peakTFMotifList <- lapply(peakTFMotifList,function(d) {d$FDR <- p.adjust(d$p.value,method="fdr"); d})

# Find diff peaks using cut-off
numsig <- unlist(lapply(peakTFMotifList,function(x) sum(x$FDR < FDR,na.rm = TRUE))) # Won't consider NAs
numsig

# Plot number of sig peaks per motif
motif.d <- data.frame("Motif"=extractTFNames(names(peakTFMotifList)),"NumSigPeaks"=numsig)
motif.d <- motif.d %>% arrange(desc(NumSigPeaks))
motif.d$Motif <- factor(as.character(motif.d$Motif),levels=as.character(motif.d$Motif))

gNumDiffPeaks <- ggplot(motif.d,aes(x=Motif,y=NumSigPeaks)) +
  geom_point(stat="identity",fill="slategray",shape=21,color="black") + theme_bw() +
  theme_classic()+
  theme(axis.text.x = element_text(angle=90,hjust=1,vjust=0.5,size=5),axis.text = element_text(color="black"))+
  labs(y=paste0("No. of differential peaks\n (FDR < ",FDR,")"))+
  scale_y_continuous(limits=c(0,max(motif.d$NumSigPeaks)+1000))
gNumDiffPeaks


# Get union of all sig peaks in sig motifs
sigPeaksIndicesPerMotif <- lapply(peakTFMotifList,function(l,FDR.cut=FDR) { which(l$FDR < FDR.cut)}) # Won't consider NAs

names(sigPeaksIndicesPerMotif)[1]

# Pooling sig peaks
sigPeaks <- Reduce("union",sigPeaksIndicesPerMotif) # Note here that taking union doesn't automatically sort it numerically (e.g. peak 2 before peak 3; we do this only in the end to save it separately)

length(sigPeaks)
sum(is.na(sigPeaks))

Zgroups <- .binarizeZMat(Z)

LFC.mat <- Matrix(0,nrow = length(sigPeaks),ncol = nrow(Z))
colnames(LFC.mat) <- extractTFNames(rownames(Z))

# Mean-normalize counts before taking fold-change
SE.norm <- centerCounts(SE)

# Takes a while to run, make sure we save the output for future ref
# Loop through each motif
for(m in 1:nrow(Zgroups)){
  # Mean accessibility for motif high cells
  countsHigh <- Matrix::rowMeans(assay(SE.norm)[sigPeaks,Zgroups[m,]==1] +1) # Add pseudocount
  # Mean accessibility for motif low cells
  countsLow <- Matrix::rowMeans(assay(SE.norm)[sigPeaks,Zgroups[m,]==0] +1) # Add pseudocount
  LFC <- log2(countsHigh/countsLow)

  LFC.mat[,m] <- LFC
}

sum(is.na(LFC.mat))

# Save to avoid re-run
saveRDS(LFC.mat,"./diffPeaks/LFC_mat_sigPeaks_jackSigMotifs.rds")

# Cluster sig peaks into modules using KNN + Louvain of log-fold change
set.seed(123)
knn <- FNN::get.knn(t(scale(t(as.matrix(LFC.mat)))), algo="kd_tree", k =30)[["nn.index"]]
igraphObj <- igraph::graph_from_adjacency_matrix(igraph::get.adjacency(igraph::graph.edgelist(data.matrix(reshape2::melt(knn)[,c("Var1", "value")]), directed=FALSE)), mode = "undirected")

# Louvain
clusters <- igraph::cluster_louvain(igraphObj)
Kmemberships <- igraph::membership(clusters)

table(Kmemberships)

# Plot number of peaks per module
K <- max(Kmemberships)
sigPeaks.d <- data.frame("Module"=1:K,"numSigGenes"=as.numeric(table(Kmemberships)))

gNumModulePeaks <- ggplot(sigPeaks.d,aes(x=factor(Module),y=numSigGenes)) + geom_bar(stat="identity",color="white",size=0.5) + theme_bw() +
  theme(axis.text.x = element_text(angle=90,hjust=1,vjust=0.5),axis.text = element_text(color="black"),axis.text.y = element_text(size=6))+
  labs(y=paste0("No. of differential peaks\n (FDR < ",FDR,")"),x="Module")+
  scale_y_continuous(expand=c(0,0),limits=c(0,max(sigPeaks.d$numSigGenes) + 500))
gNumModulePeaks
ggsave(filename = "./figs/ModulesnumDiffPeaks.pdf",height=4,width=4)


# For reference, save the derived K module maps for each reference peak (vector)
names(Kmemberships) <- paste0("Peak ",sigPeaks)

# Now we can cluster each K group separately
hlist <- list()
# Instead of 1-11, we can use a custom order after later inspection
clusterGroupOrder <- 1:K
for(i in 1:max(K)){
  k <- clusterGroupOrder[i]
  cat("Clustering peaks from K=",k," ..\n")
  kmat <- as.matrix(LFC.mat[Kmemberships %in% k,])
  # Correlation based clustering
  d <- as.dist(1-cor(t(kmat), method="pearson"))
  h <- hclust(d)
  # Return matrix of values row-sorted by clustering order (no column re-ordering)
  hlist[[i]] <- data.frame(kmat[h$order,])
}

# Merge K-module ordered, individually clustered list of matrices into single matrix
p <- dplyr::bind_rows(hlist)

# Define memership splits for K-ordered heatmap to split heatmap rows by K module
Kmemberships.sorted <- factor(rep(clusterGroupOrder,sapply(hlist,nrow)),levels = clusterGroupOrder) # Levels define the order in which split heatmap groups are plotted

# Plot heatmap based on specified ordering, split and clustered per module
library(ComplexHeatmap)
heat <- Heatmap(t(scale(t(p))),name = "log2 FC \naccessibility",
                cluster_rows = FALSE,
                column_names_gp = gpar(fontsize = 6.5),
                clustering_distance_columns = "pearson",
                col = jdb_palette("brewer_jamaica"),
                split = Kmemberships.sorted,
                show_row_names = FALSE,
                use_raster = TRUE)

draw(heat)



# Build annotation matrix of peaks x K modules
# This annotation matrix is what is used as input to chromVAR to score single cells for modules
Kannot <- sparseMatrix(i = sigPeaks, # Indices for the significant peaks returned
                       j = Kmemberships, # K cluster assignments for each of the same peaks
                       dims=c(nrow(SE),K)) # Dimensions pertaining to final matrix

colnames(Kannot) <- paste0("K",1:K)