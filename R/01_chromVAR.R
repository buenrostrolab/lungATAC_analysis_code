### Author: Vinay Kartha
### Contact: <vinay_kartha@g.harvard.edu>
### Affiliation: Buenrostro Lab, Department of Stem Cell and Regerative Biology, Harvard University

library(chromVAR)
library(chromVARmotifs)
library(motifmatchr)
library(BSgenome.Mmusculus.UCSC.mm10)
library(SummarizedExperiment)
library(tictoc)
BiocParallel::register(BiocParallel::MulticoreParam(4, progressbar = TRUE))

setwd("<data_analysis_folder>")

# Load filtered scATAC peak counts object
SE <- readRDS("./atac.se.rds")

# Get GC content for chromVAR
SE <- addGCBias(SE,genome=BSgenome.Mmusculus.UCSC.mm10)

# Peaks with zero accessibility across cells
table(Matrix::rowSums(assay(SE))==0)

# Get background peaks for chromVAR
bg <- getBackgroundPeaks(object = SE,niterations=250)
saveRDS(bg,"./chromVAR/bg_peaks.rds")

# K-mer scoring (used for clustering single cells)
kmer_ix <- matchKmers(k = 6,subject = SE,genome=BSgenome.Mmusculus.UCSC.mm10)
saveRDS(kmer_ix,"./chromVAR/kmer_ix.rds")

# Motif scoring (used to annotate single cells)
motif_ix <- matchMotifs(pwms = mouse_pwms_v2,subject = SE,genome=BSgenome.Mmusculus.UCSC.mm10)
saveRDS(motif_ix,"./chromVAR/motif_ix.rds")

tic("Computing deviations for motifs ..\n\n")
dev_motif <- computeDeviations(object = SE,annotations = motif_ix,background_peaks=bg)
saveRDS(dev_motif,"./chromVAR/devMotif.rds")
toc()

tic("Computing deviations for k-mers ..\n\n")
dev_kmer <- computeDeviations(object = SE,annotations = kmer_ix,background_peaks=bg)
saveRDS(dev_kmer,"./chromVAR/devKmer.rds")
toc()
