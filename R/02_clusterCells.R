library(SummarizedExperiment)
library(chromVAR)
library(BuenRTools)
library(BuenColors)
library(ggplot2)
library(leiden)
library(densityClust)

splitAndFetch <- function(vec,
                          delim,
                          part){
  if(length(part)==1){
    sapply(strsplit(as.character(vec),delim),"[[",part) } else {
      sapply(strsplit(as.character(vec),delim),function(x) paste(x[part],collapse = delim))
    }
}

myShuf <- function(df,seed=123){
  set.seed(seed)
  BuenColors::shuf(df)
}


################################## LOAD DATA #################################################

setwd("/mnt/users/vinay/analysis/lungATAC/results/manuscript_dataFreeze/resubmission/")

# Load merged SE object (contains final filtered cells: see /mnt/users/vinay/analysis/lungATAC/results/manuscript_dataFreeze/resubmission/code/filtandMergeSE.R)
SE <- readRDS("./SE/SE_merged_filt.rds")
table(SE$TypeNew)

dev.kmer <- readRDS("./chromVAR/devKmer.rds")
dev.motif <- readRDS("./chromVAR/devMotif.rds")
# These are the original submission module scores
dev.mod <- readRDS("./chromVAR/devK10_diffpeaks.rds")

stopifnot(all.equal(colnames(SE),colnames(dev.kmer)) & all.equal(colnames(SE),colnames(dev.motif)) & all.equal(colnames(SE),colnames(dev.mod)))

# Some of these motifs are redundant (since it was made for mouse_pwms_v2)
# Filter in pwm_v3 motifs (797 of them)
pwmv3 <- load("/mnt/users/vinay/data/annot/motifs/mouse_pwms_v3.RData")
dev.motif <- dev.motif[names(mouse_pwms_v3),]
Zmotif <- deviationScores(dev.motif)
which(is.na(Zmotif),arr.ind = TRUE)

# Set these to 0, and save mat object
Zmotif[is.na(Zmotif)] <- 0
R.matlab::writeMat(Z=Zmotif,"./forJason/Zmotif.mat") # Contains all cells, all motifs

# We could set these to 0, but instead, let us bag the motifs, since those don't have 0s in them

dev.motif.bagged <- bagDeviations(object = dev.motif,organism = "mouse",cor = 0.7)
which(is.na(deviationScores(dev.motif.bagged)),arr.ind = TRUE)
saveRDS(dev.motif.bagged,"./chromVAR/devMotif_bagLeaders.rds")

Zmotif <- deviationScores(dev.motif.bagged)
sum(is.na(Zmotif)) # Check
rownames(Zmotif) <- extractTFNames(rownames(Zmotif))

# Export bagged motif matrix for Jason as .mat file
R.matlab::writeMat(Z=Zmotif,"./forJason/Zmotifbagged.mat") # Contains all cells, bagged motifs

Zk <- deviationScores(dev.kmer)
sum(is.na(Zk))

###############################################################################################
############################# TUMOR MET AND NORMAL DATA CLUSTERING ############################

# Run PCA on all data (exclude ETP cells here)
ETPcells <- grepl("Early",SE$TypeNew,fixed = FALSE)
table(SE$TypeNew,ETPcells) # Check

# If running PCA LEAVING OUT ETP CELLS
pc <- cachePCA(dataSet=t(scale(Zk[,-which(ETPcells)])),
               scale=FALSE,
               center=FALSE,
               cachePath = "./PCAcache/")

pc.scores <- pc$x[,1:20] # Top20 PCs
saveRDS(pc.scores,"./chromVAR/k6mer_PCscores_noETP.rds")

# This uses  the UWOT packagewhich
# umap <- cacheUMAP2(dataSet = pc.scores,
#                    cachePath = "./UMAPcache/",
#                    metric="cosine",
#                    seed=200,
#                    useUWOT = TRUE,
#                    n_neighbors=20,
#                    min_dist=0.5,
#                    ret_model=TRUE,
#                    runAnyway = TRUE) # Have to set this with UMAP, re-loading it won't let you do projection, if needed

# Changed parameters
umap <- cacheUMAP2(dataSet = pc.scores,
                   cachePath = "./UMAPcache/",
                   metric="cosine",
                   seed=123,
                   useUWOT = TRUE,
                   n_neighbors=20,
                   min_dist=0.4,
                   ret_model=TRUE,
                   runAnyway = TRUE) 
saveRDS(umap,"./chromVAR/k6mer_UMAPobj_noETP.rds")

umap.d <- as.data.frame(umap$embedding)
colnames(umap.d) <- paste0("UMAP",c(1,2))
rownames(umap.d) <- rownames(pc.scores)

umap.d$Exp <- SE$TypeNew[-which(ETPcells)]
umap.d$TissueType <- ifelse(grepl("Met",umap.d$Exp),"Met",
                            ifelse(grepl("Normal",umap.d$Exp),"Normal",
                                   ifelse(grepl("Early",umap.d$Exp),"Early tumor",
                                   "Tumor")))
# Base
g1 <- ggplot(umap.d,aes(x=UMAP1,y=UMAP2)) + 
  geom_point(size=0.1,color="darkslategray") + 
  theme_classic()

g1

# With tissue type
g2 <- ggplot(umap.d,aes(x=UMAP1,y=UMAP2,color=TissueType)) + 
  geom_point(size=0.1) + 
  theme_classic() + 
  scale_color_manual(values=c("firebrick","steelblue","gray")) + 
  guides(colour = guide_legend(override.aes = list(size=5)))

g2

# Flip y-axis (and reverse umap1)
umap.d$UMAP1 <- -umap.d$UMAP1

umap.d$TissueType <- as.character(umap.d$TissueType)
g2.flipped <- ggplot(myShuf(umap.d,seed = 1),aes(x=UMAP1,y=UMAP2,color=TissueType)) + 
  geom_point(size=0.01,alpha=0.9) +
  theme_classic() + 
  scale_color_manual(values=c("Tumor"="gray","Met"="firebrick","Normal"="steelblue")) + 
  theme(legend.position="none",legend.key.size = unit(0, 'lines'),
        legend.text = element_text(size=5),legend.title = element_text(size=8),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),axis.title = element_blank())+
  guides(colour = guide_legend(override.aes = list(size=1))) + 
  coord_flip() # NOTE WE'RE PLOTTING UMAP2 VS UMAP1!

g2.flipped

ggsave(filename = "./figures/UMAP_tissueType_rotated.png",plot = g2.flipped,width = 4.5,height = 3.5)
ggsave(filename = "./figures/UMAP_tissueType_rotated.pdf",plot = g2.flipped,width = 4.5,height = 3.5)

# Color by FRIP/depth
umap.d$FRIP <- SE$FRIP[-which(ETPcells)]
umap.d$depth <- log10(SE$depth[-which(ETPcells)])
library(ggrastr)
g4a <- ggplot(umap.d,aes(UMAP1,UMAP2,color=FRIP)) + 
  geom_point_rast(size=0.01) + 
  scale_color_gradientn(colours = jdb_palette("flame_light")) +  
  theme_classic() +
  theme(legend.text = element_text(size=4),
        legend.title = element_text(size=7),
        legend.key.size = unit(0.15,"in"),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),axis.title = element_blank())+
  coord_flip()


g4b <- ggplot(umap.d,aes(UMAP1,UMAP2,color=depth)) + 
  geom_point_rast(size=0.01) + 
  scale_color_gradientn(colours = jdb_palette("flame_light")) +  
  theme_classic() +
  theme(legend.text = element_text(size=4),
        legend.title = element_text(size=7),
        legend.key.size = unit(0.15,"in"),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),axis.title = element_blank())+
  coord_flip()
g4comb <- cowplot::plot_grid(g4a,g4b,nrow=1)

ggsave(filename = "./figures/UMAP_QC_rotated.pdf",plot = g4comb,width = 8,height = 3.5)
ggsave(filename = "./figures/UMAP_QC_rotated.png",plot = g4comb,width = 8,height = 3.5)



# KNN based on pc scores (used for gene score smoothing etc.)
set.seed(123)
knn <- FNN::get.knn(pc.scores, algo="kd_tree", k = 20)[["nn.index"]]
rownames(knn) <- rownames(pc.scores)
saveRDS(knn,"./chromVAR/k6mer_PCscores_noETP_kNN20.rds")

# Plot smoothed gene scores for select markers
# Load raw TSS counts for all 
TSS <- readRDS("./TSSscores/TSS_merged_filt_norm.rds")
all.equal(colnames(SE),colnames(TSS))
markerlist <- c("Cebpa","Sftpb","Sftpc","Bmp1","Bmp2","Cav1","Fgf18","Mir7680","Timp3","Ptprc","Tead4","Nkx2-4","Runx2","Fosl1","Igfbp3","Pdgfra")
all(markerlist %in% rownames(TSS))

# Smooth TSS scores for these markers
library(doParallel)
TSSsmoothed <- smoothGeneScoresNN3(NNmat = knn,TSSmat = TSS[,-which(ETPcells)],geneList = markerlist,nCores = 6)

gglist <- list()
for(gene in markerlist){
gglist[[gene]] <-   plotMarker2D(df = umap.d[,2:1], # UMAP2 VS UMAP1
                              markerScore = TSSsmoothed[gene,],
                              pointSize = 0.1,
                              maxCutoff = "q0.99",
                              rasteRize = FALSE,
                              colorPalette = "wolfgang_basic",
                              markerName = "Genescore",
                              plotTitle = gene)
}

cowplot::plot_grid(plotlist=gglist,nrow=3)
###############################################################################################

############################ NORMAL CELL CLUSTERING AND ANNOTATION ############################
normalCells <- grepl("Normal",SE$TypeNew,fixed = FALSE)
table(SE$TypeNew,normalCells)
# RUN PCA ON K-MER SCORES USING THE NORMALS ONLY
pcNormal <- cachePCA(dataSet=t(scale(Zk[,normalCells])),
                     scale=FALSE,
                     center=FALSE,
                     cachePath = "./PCAcache/")

pc.scoresNormal <- pcNormal$x[,1:20] # Top20 PCs
saveRDS(pc.scoresNormal,"./chromVAR/k6mer_PCscores_normalCells.rds")

# Louvain clustering
set.seed(123)
knnNormal <- FNN::get.knn(pc.scoresNormal, algo="kd_tree", k = 50)[["nn.index"]]
rownames(knnNormal) <- colnames(SE)[normalCells]
saveRDS(knnNormal,"./chromVAR/k6mer_PCscores_normalCells_kNN50.rds")

igraphObj <- igraph::graph_from_adjacency_matrix(igraph::get.adjacency(igraph::graph.edgelist(data.matrix(reshape2::melt(knnNormal)[,c("Var1", "value")]), directed=FALSE)), mode = "undirected")
KmembershipsLouvain <- igraph::membership(igraph::cluster_louvain(igraphObj))
table(KmembershipsLouvain)

# UMAP to visualize
# This uses  the UWOT packagewhich
umapNormal <- cacheUMAP2(dataSet = pc.scoresNormal,
                   cachePath = "./UMAPcache/",
                   metric="cosine",
                   n_neighbors=20,
                   min_dist=0.5)

umapNormal.d <- as.data.frame(umapNormal)
colnames(umapNormal.d) <- paste0("UMAP",c(1,2))
rownames(umapNormal.d) <- rownames(pc.scoresNormal)

normalCellBarcodes <- splitAndFetch(rownames(umapNormal.d),"_",2) # This is to match the barcode style used originally (without experiment added)
umapNormal.d$Louvain <- KmembershipsLouvain

gNormLouvain <- ggplot(umapNormal.d,aes(x=UMAP1,y=UMAP2,color=factor(Louvain))) + 
  geom_point(size=0.1) + 
  theme_classic() + 
  scale_color_manual(values=jdb_palette("lawhoops"))+
  guides(colour = guide_legend(override.aes = list(size=3)))+
  labs(color="Louvain")
gNormLouvain

# Louvain looks better so we keep it for now
umapNormal.d$myLabel <- paste("Normal",umapNormal.d$Louvain)
clustcentroids <- umapNormal.d %>% group_by(Louvain,myLabel) %>% summarise(UMAP1=median(UMAP1),UMAP2=median(UMAP2))
rownames(clustcentroids) <- clustcentroids$myLabel
dend <- hclust(dist(as.matrix(clustcentroids[,3:4])))
plot(dend)

library(ggrepel)
my.cols <- colorRampPalette(c("firebrick3","orange","palegreen4","lightskyblue1","navy"))(max(KmembershipsLouvain))
names(my.cols) <- dend$order

gNormLouvainID <- gNormLouvain + theme(legend.position="none")+
  scale_color_manual(values = my.cols)+
  labs(color="Louvain")+
  geom_label_repel(data=clustcentroids,
                   aes(x=UMAP1,y=UMAP2,fill=factor(Louvain),label=myLabel),
                   color="white",show.legend = FALSE,size=3,alpha=0.95,segment.colour = "darkgray",segment.size = 1) + 
  scale_fill_manual(values=my.cols)

cowplot::plot_grid(gNormLouvain,gNormLouvainID,nrow=1)

# Match these to what we generated before
# This is metadata associated with the original submission
meta <- read.table("../tSNE/tSNE_coords_with_info.txt",sep="\t",stringsAsFactors = FALSE,header = TRUE)
dim(meta)
table(meta$clusterLouv)
all(normalCellBarcodes %in% rownames(meta))

# Fetch the original Louvain cluster for each of these normal barcodes
umapNormal.d$Louvainoriginal <- meta[normalCellBarcodes,"clusterLouv"]

x <- table(umapNormal.d$Louvain,umapNormal.d$Louvainoriginal)
# Here, we look for the max overlap in cluster associations (i.e. determining which old cluster corresponds to which new cluster)
analog <- apply(x,1,function(a) { names(which.max(a))})
names(analog) <- paste("Normal",names(analog))

# Load the annotation label
normAnnot <- read.csv("./annot/normalLouvClusterAnnotations_Lafave.csv",header=FALSE,row.names = 1,col.names = c("Cluster","Label"),stringsAsFactors = FALSE)
newlabs <- normAnnot[analog,]
names(newlabs) <- names(analog)

# Add annotation to cluster data
umapNormal.d$Label <- newlabs[umapNormal.d$Louvain]
clustcentroids$Label <- newlabs

# Update colors
normcols <- readRDS("../tSNE/clusterColPalette_NEW.rds")
newCols <- normcols[analog]
names(newCols) <- 1:length(newCols)

# Modify color for Monocytes and Macrophages II/III cluster
# This is Louvain cluster # 1 and 2
newCols[1:2] <- c("#387A63","#6EC5A7")

gNormLouvainAnnot <- gNormLouvain + theme(legend.position="none")+
  scale_color_manual(values = newCols)+
  labs(color="Louvain")+
  geom_label_repel(data=clustcentroids,
                   aes(x=UMAP1,y=UMAP2,fill=factor(Louvain),label=Label),
                   color="white",show.legend = FALSE,size=3,alpha=0.95,segment.colour = "darkgray",segment.size = 1) + 
  scale_fill_manual(values=newCols)

gNormLouvainAnnot

ggsave(plot = gNormLouvainAnnot,filename = "./figures/UMAPnormals.pdf",width = 6,height = 5,useDingbats=FALSE)

# Make meta data table with cluster label annotations for normals, in addition to tumor and mets
umap.d$Louvain <- umap.d$TissueType
umap.d$Annotation <- umap.d$TissueType

# Update normal annotation Louvain cluster #s
umap.d[rownames(umapNormal.d),"Louvain"] <- umapNormal.d$myLabel # Labels like "Normal 1"
umap.d[rownames(umapNormal.d),"Annotation"] <- umapNormal.d$Label # Labels like "Monocytes"

# Just saving two copies
write.table(umap.d,"./forJason/cellMetaNoETP.txt",sep="\t",row.names = TRUE,col.names = TRUE,quote = FALSE)
write.table(umap.d,"./annot/cellMetaNoETP.txt",sep="\t",row.names = TRUE,col.names = TRUE,quote = FALSE)

# Also write one with just normal barcodes and their annotation (to make split bams for tracks)
write.table(data.frame("Barcode"=splitAndFetch(rownames(umap.d)[umap.d$TissueType %in% "Normal"],"_",2),
                       "Cluster"=gsub(x=gsub(x=umap.d[umap.d$TissueType %in% "Normal","Annotation"],pattern = " ",replacement = ""),pattern = "/",replacement = "_")),
                  file="./forJason/normalBarcodesLouvainClust.txt",sep="\t",row.names = FALSE,col.names = FALSE,quote=FALSE)


# Create new color scheme, including tumor and mets, and ETPs
majorCols <- newCols
names(majorCols) <- newlabs
majorCols <- c(majorCols,"Tumor"="lightgray","Met"="#B02725")
saveRDS(majorCols,"./annot/myCols.rds")

all(markersToPlot %in% rownames(TSSsmoothed))

TSSsmoothedNorm <- smoothGeneScoresNN3(NNmat = knnNormal,
                                   TSSmat = TSS[,normalCells],
                                   geneList = NULL,
                                   nCores = 6)


all.equal(rownames(umapNormal.d),colnames(TSSsmoothedNorm))

# For marker genes, create Seurat object and find top markers unique to each Louvain cluster
library(Seurat)
normGS <- CreateSeuratObject(counts = TSSsmoothedNorm,assay = "ATAC")
normGS$Louvain <- umapNormal.d$Label
normGS <- SetIdent(object = normGS,value = "Louvain")
normMacMarkers <- FindMarkers(normGS,ident.1 = "Macrophages I",only.pos = TRUE)

markerList <- list()
pdf("./figures/UMAPnormalsGeneScores.pdf",width = 5,height = 5,useDingbats = FALSE)
for(gene in markersToPlot){
markerList[[gene]] <- plotMarker2D(df = umapNormal.d[,1:2],
                                   markerScore = TSSsmoothed[gene,],plotTitle = gene,
                                   markerName = "Genescore",
                                   rasteRize = FALSE,
                                   maxCutoff = "q0.99",
                                   colorPalette = "brewer_fire",
                                   pointSize = 0.1)
print(markerList[[gene]])
}

dev.off()


################################################################################

#################### TUMOR + MET CELL FILTERING ###############################
# Perform density clustering in order to remove outlier cells
# Density clustering of all data
umapdist = dist(umap.d[,1:2])
set.seed(0)
dclust = densityClust(umapdist,gaussian=T)

options(repr.plot.width=6, repr.plot.height=6)

plot(dclust$rho,dclust$delta,pch=20,cex=0.6)


# IMPORTANT NOTE: SINCE WE CHANGED THE UMAP LATER, THIS WILL AFFECT THE DENSITY CLUSTERING MODERATELY
# APPLYING THESE PARAMETERS TO THE CHANGED UMAP YIELDED 6 EXTRA CELLS IN THE MET+TUMOR FILTERED BARCODE SPACE OBTAINED BEFORE
# WE DID NOT RE-RUN MODULES SINCE THIS WOULD MAKE A MINOR DIFFERENCE.
# NOTE FOR RECORD KEEPING THAT THE TUM+MET CELL BARCODES SAVED WERE NOT OVERWRITTEN WITH THE NEW DENSITY CLUSTERING FILTER
rho_cutoff <- 50
delta_cutoff <- 0.5

dclust = findClusters(dclust, rho = rho_cutoff, delta = delta_cutoff)

plot(dclust$rho,dclust$delta,pch=20,cex=0.6)
points(dclust$rho[dclust$peaks],dclust$delta[dclust$peaks],col="red",pch=20,cex=0.8)
text(dclust$rho[dclust$peaks]-2,dclust$delta[dclust$peaks]+0.5,labels=dclust$clusters[dclust$peaks])


umap.d$clusterDens <- dclust$clusters

gDens <- ggplot(umap.d,aes(x=UMAP1,y=UMAP2,color=factor(clusterDens))) +
  geom_point(size=0.2) + scale_color_manual(values=jdb_palette("lawhoops")) + theme_classic() + coord_flip()

# Add centroid labels to UMAP based on density clustering
centroids <- umap.d %>% group_by(clusterDens) %>% summarise(UMAP1=median(UMAP1),UMAP2=median(UMAP2))

library(ggrepel)
gDensLab <- ggplot(umap.d,aes(x=UMAP1,y=UMAP2,color=factor(clusterDens))) + 
  geom_point(size=0.1) + theme_classic() + 
  scale_color_manual(values=jdb_palette("lawhoops")) + 
  geom_label_repel(data=centroids,aes(x=UMAP1,y=UMAP2,label=clusterDens,fill=factor(clusterDens)),
                   color="white",show.legend = FALSE,size=2)+
  scale_fill_manual(values=jdb_palette("lawhoops"))+
  theme(legend.position = "none") + coord_flip()

ggsave(plot = gDensLab,filename = "./figures/UMAP_densityClust_rotated.pdf",height = 4,width = 4.5)

# Plot density cluster by experiment stacked bar plot (fraction of cells)
m <- as.matrix(table(umap.d$TissueType,umap.d$clusterDens))

# Convert to fraction
m <- t(t(m)/colSums(m))

library(ComplexHeatmap)
pdf("./figures/Heatmap_densityClust_tissueType.pdf",width = 5,height = 2)
Heatmap(as.data.frame.matrix(m),name = "Fraction\nof cells",
        col = jdb_palette("flame_light"),
        row_names_gp = gpar(fontsize = 10),
        column_names_gp = gpar(col=jdb_palette("lawhoops")),
        clustering_distance_columns = "pearson",
        clustering_distance_rows = "pearson")
dev.off()

# Cells to Keep for module analysis (and zoomed in tumor/met plots)
# Define by both density clustering and experiment (to remove sorting outliers)
clustersToRemove <- c(3:6,9:10,13,14,17,18) # Note these clusters will not cover AT1/AT2 since those overlap tumor cells
# This tumor and met filter will now filter out AT1/AT2 cells
cellsToKeep <- rownames(umap.d)[umap.d$TissueType %in% c("Tumor","Met")  & !(umap.d$clusterDens  %in% clustersToRemove)]

# Check cells chosen for modules
gTumMet <- ggplot(umap.d,aes(UMAP1,UMAP2)) + geom_point(size=0.1,color="lightgray") + theme_classic() + 
  geom_point(data=umap.d[cellsToKeep,],color="darkorange",size=0.3) + coord_flip()

#saveRDS(cellsToKeep,"./modules/tumorMetCellsModules.rds")

ggsave(plot = gTumMet,filename = "./figures/UMAP_MetTumChosen_highlight_rotated.pdf",height = 4,width = 4.5)

# Plot UMAP zoomed in with only tumor / met cells
gZoom <- ggplot(umap.d[cellsToKeep,],aes(x=UMAP1,y=UMAP2,color=Annotation))+
  geom_point(size=0.01)+
  theme_classic()+
  scale_color_manual(values = majorCols) + 
  theme(legend.position="none",legend.key.size = unit(0, 'lines'),
        legend.text = element_text(size=5),legend.title = element_text(size=8),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),axis.title = element_blank())+
  guides(colour = guide_legend(override.aes = list(size=1))) + 
  coord_flip()

gZoom

# Plot for specific tumors
# Also do it for the remaining tumors (supplement)
tumorsToPlot <- c("Tumor-M2T3","Tumor-M3T2","Tumor-M1T1", "Tumor-M1T2" ,"Tumor-M1T3", "Tumor-M3T1",  "Tumor-M3T3")
all(tumorsToPlot %in% umap.d$Exp)

gTlist <- list()
pdf("./figures/UMAP_TumorsHighlight_rotated.pdf",height = 4,width=4.5)
for(i in tumorsToPlot){
  
  # Highlight certain tumors
  gTlist[[i]] <- ggplot(umap.d,aes(x=UMAP1,y=UMAP2)) + 
    geom_point(size=0.01,alpha=0.9,color="gray") +
    geom_point(data=umap.d[umap.d$Exp %in% i,],size=0.01,color="firebrick")+
    labs(title=i)+
    theme_classic() + 
    theme(legend.position="none",legend.key.size = unit(0, 'lines'),
          legend.text = element_text(size=5),legend.title = element_text(size=8),
          axis.line = element_blank(),
          axis.ticks = element_blank(),
          axis.text = element_blank(),
          axis.title = element_blank(),plot.title = element_text(hjust=0.5))+
    coord_flip() # NOTE WE'RE PLOTTING UMAP2 VS UMAP1!
  print(gTlist[[i]])
}
dev.off()


##############################################################################

#################### EARLY TIME POINT PROJECTION #############################

#Predict ETP coords using embedding of rest of data
library(uwot)
# First let us get projected PC scores using k-mer Z scores for the ETP cells
# Here, we use the PCA run for single cells on k-mers where we left out the ETP
ETPPCscores <- t(Zk[,ETPcells]) %*% pc$rotation

# Now, we can get umap coords for the top 20 PCs (same used for non-ETP cells) for these ETPs
ETPcoords <- as.data.frame(umap_transform(ETPPCscores[,1:20],umap))

rownames(ETPcoords) <- colnames(SE[,ETPcells])
colnames(ETPcoords) <- c("UMAP1","UMAP2")
ETPcoords$UMAP1 <- -(ETPcoords$UMAP1) # This is because we also reverse the original UMAP1 coords

ggplot(umap.d,aes(x=UMAP1,y=UMAP2)) + 
  geom_point(size=0.1,color="darkslategray") + 
  theme_classic() + geom_point(data = ETPcoords,size=0.05,color="gold") + coord_flip() # Flip coords as before (UMAP2 vs UMAP1)

# Save ETP coords
write.table(ETPcoords,"./forJason/ETP_projected_UMAP_coords.txt",sep="\t",quote=FALSE)
