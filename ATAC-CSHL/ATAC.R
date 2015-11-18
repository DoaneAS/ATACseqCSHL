library(Rsamtools)
library(Biostrings)
# install using biocLite if needed

# set wd as needed
setwd("TEACHING")

# set bamfile variable
bamfile <- 'LY1_ATAC_chr18.bam'

# scan it into mem
bam <- scanBam(bamfile)

# get fragment size directly from BAM file
fragsize <- bam[[1]]$isize [ bam[[1]]$isize > 1 & bam[[1]]$isize < 500]

# look at a few values
head(fragsize)

# plot fragment size density
plot(table(fragsize))

library(ChIPpeakAnno)

# load up BED file with ATAC peaks (MACS2 derived)
peaks <- toGRanges("ATAC_LY1_peaks_chr18.bed")

# EXERCISE:how many peaks ? average peak width ?

# next steps: annotate peaks with genes

# attach TSS data from human genes (GRCh37/hg19)
data(TSS.human.GRCh37)

# load up package
# biocLite("ChIPpeakAnno")
library(ChIPpeakAnno)

# annotate peaks with TSS (window of 2kb)
anno <- annotatePeakInBatch(peaks, AnnotationData=TSS.human.GRCh37,
                            output="overlapping", maxgap=1000L)

# look at it
anno

# load gene name database
biocLite("org.Hs.eg.db")
library(org.Hs.eg.db)

# add gene names
anno <- addGeneIDs(annotatedPeak=anno,
                        orgAnn="org.Hs.eg.db",
                        IDs2Add="symbol")

anno

# count 
table(as.data.frame(anno)$insideFeature)

# make pie chart
pie(table(as.data.frame(anno)$insideFeature))

# load up PU.1 ChIPseq peaks
PU1peaks <- toGRanges("PU1_LY1_peaks_chr18.bed")

# how many peaks
length(PU1peaks)

# annotate ATAC peaks with PU.1 peaks
annoPU1 <- annotatePeakInBatch(peaks, AnnotationData=PU1peaks, output="overlapping")

# another way 
ol <- findOverlapsOfPeaks(peaks, PU1peaks, maxgap=0)

# get lists of peaks
peaklist <-ol$peaklist

names(peaklist)

# different classes of overlap/nonoverlap
length(peaklist[["peaks///PU1peaks"]])
length(peaklist[["peaks"]])
length(peaklist[["PU1peaks"]])

# visual
makeVennDiagram(ol, totalTest=1e+5)

# get sequence data for hg19
biocLite("BSgenome.Hsapiens.UCSC.hg19")
library("BSgenome.Hsapiens.UCSC.hg19")

# get peaks
peaksWithSequences <- getAllPeakSequence(peaks, upstream=0,
                                         downstream=0, genome=Hsapiens)

# look
peaksWithSequences

# save as DNAStringSet
s <- DNAStringSet(peaksWithSequences$sequence)

# load up motif database
biocLite("MotifDb")
biocLite("seqLogo")
library(seqLogo)
library(MotifDb)

# get PU1 motif 
pfm.pu1 <- query(MotifDb,"Spib")[[1]]

# look at it
pfm.pu1

# get seq logo (needs seqLogo package)
seqLogo(pfm.pu1)

# to integer
pcm.pu1 <- round(100 * pfm.pu1)

# motif matching
mmpu1 <- matchPWM(pcm.pu1, unlist(s), "90%")

#
length(mmpu1)

# get peak sequence for shifted/adjacent
peaks_shifted = shift(peaks, width(peaks))
peaksWithSequences_shifted <- getAllPeakSequence(peaks_shifted, upstream=0,
                                         downstream=0, genome=Hsapiens)

# get DNAStringSet
ss <- DNAStringSet(peaksWithSequences_shifted$sequence)

# match to adjacent seq
mmpu1_s <- matchPWM(pcm.pu1, unlist(ss), "90%")

# how many motifs
length(mmpu1_s)

# EXERCISE: what is the significance of the true number of motif matches ?

# EXERCISE: how many reads under the ATAC peaks



