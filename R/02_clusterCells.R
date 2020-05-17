### Author: Vinay Kartha
### Contact: <vinay_kartha@g.harvard.edu>
### Affiliation: Buenrostro Lab, Department of Stem Cell and Regerative Biology, Harvard University


library(SummarizedExperiment)
library(chromVAR)
library(uwot)
library(ggplot2)


setwd("<data_analysis_folder>")

# Load filtered scATAC peak counts object
SE <- readRDS("./atac.se.rds")

# Load k-mer deviation scores (see 01_chromVAR.R)
dev.kmer <- readRDS("./chromVAR/devKmer.rds")
Zk <- deviationScores(dev.kmer)

# For the initial clustering, we leave out the Early Time Point (ETP) cells
# Run PCA on this subset first

ETPcells <- grepl("Early",SE$TypeNew,fixed = FALSE)
table(SE$TypeNew,ETPcells) # Check

# Run PCA using k-mer Z scores
pc <- prcomp(x=t(scale(Zk[,-which(ETPcells)])),
               scale=FALSE,
               center=FALSE)

pc.scores <- pc$x[,1:20] # Top20 PCs

# Run UMAP projection
set.seed(seed=123)
umap <- umap(X = pc.scores,
             metric="cosine",
             n_neighbors=20,
             min_dist=0.4,
             ret_model=TRUE)


umap.d <- as.data.frame(umap$embedding)
colnames(umap.d) <- paste0("UMAP",c(1,2))
rownames(umap.d) <- rownames(pc.scores)

# Meta-data
umap.d$Exp <- SE$TypeNew[-which(ETPcells)]
umap.d$TissueType <- ifelse(grepl("Met",umap.d$Exp),"Met",
                            ifelse(grepl("Normal",umap.d$Exp),"Normal",
                                   ifelse(grepl("Early",umap.d$Exp),"Early tumor",
                                          "Tumor")))
umap.d$TissueType <- as.character(umap.d$TissueType)


# Visualize clustering

# Flip y-axis (and reverse umap1)
umap.d$UMAP1 <- -umap.d$UMAP1
ggplot(myShuf(umap.d,seed = 1),aes(x=UMAP1,y=UMAP2,color=TissueType)) +
  geom_point(size=0.01,alpha=0.9) +
  theme_classic() +
  scale_color_manual(values=c("Tumor"="gray","Met"="firebrick","Normal"="steelblue")) +
  theme(legend.position="none",legend.key.size = unit(0, 'lines'),
        legend.text = element_text(size=5),legend.title = element_text(size=8),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),axis.title = element_blank())+
  guides(colour = guide_legend(override.aes = list(size=1))) +
  coord_flip()

# Color cell by Fraction of Reads in Peaks (FRIP) and sequencing depth
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


#################################### ETP cell projection ####################################
# Multiple Z scores of ETP cells with pc embeddings from rest of data to get projected ETP cell scores
ETP_PCscores <- t(Zk[,ETPcells]) %*% pc$rotation

# Now, we can get umap coords for the top 20 PCs (same used for non-ETP cells) for these ETPs
ETPcoords <- as.data.frame(umap_transform(ETP_PCscores[,1:20],umap))

rownames(ETPcoords) <- colnames(SE[,ETPcells])
colnames(ETPcoords) <- c("UMAP1","UMAP2")
ETPcoords$UMAP1 <- -(ETPcoords$UMAP1) # This is because we also reversed the original UMAP1 coords

ggplot(umap.d,aes(x=UMAP1,y=UMAP2)) +
  geom_point(size=0.1,color="darkslategray") +
  theme_classic() + geom_point(data = ETPcoords,size=0.05,color="gold") + coord_flip() # Flip coords as before (UMAP2 vs UMAP1)

