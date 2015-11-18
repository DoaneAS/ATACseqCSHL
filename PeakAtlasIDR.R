

library(GenomicRanges)
library(GenomicAlignments)
library(rtracklayer)
source("/Volumes/AD/Projects/ATAC/bin/atacseqFuncitons.R")


extraCols_IDR <- c(peakName = "character", IDRscore = "numeric", strand="character", signalValue = "numeric", pvalue = "numeric")

#NBidr <- import.bed("NB.idr.bed", extraCols=extraCols_IDR)

#CBidr <- import.bed("CB.idr.bed", extraCols=extraCols_IDR)

NBidr <- import.bed("sample_NB_IDR.narrow.all.bed", extraCols=extraCols_IDR)

CBidr <- import.bed("sample_CB_IDR.narrow.all.bed", extraCols=extraCols_IDR)

NBidrb <- import.bed("sample_NB_IDR.broad.bed", extraCols=extraCols_IDR)

CBidrb <- import.bed("sample_CB_IDR.broad.bed", extraCols=extraCols_IDR)

seqnames(CBidr)

NBpeaks <- NBidr[NBidr$IDRscore > 400,]
CBpeaks <- CBidr[CBidr$IDRscore > 400,]

NBpeaks <- NBidr[NBidrb$IDRscore > 400,]
CBpeaks <- CBidr[CBidrb$IDRscore > 400,]


export.bed(NBpeaks, "NB.IDR.peaks.bed")
export.bed(CBpeaks, "CB.IDR.peaks.bed")


median(width(NBpeaks))

names(NBpeaks) <- paste(seqnames(NBpeaks), start(NBpeaks), end(NBpeaks), "NB", sep = ":")
NBpeaks$peakName.NB <- names(NBpeaks)

names(CBpeaks) <- paste(seqnames(CBpeaks), start(CBpeaks), end(CBpeaks), "CB", sep= ":")
CBpeaks$peakName.CB <- names(CBpeaks)





nbol.any <- subsetByOverlaps(NBpeaks, CBpeaks, minoverlap=1, ignore.strand=T) 
cbol.any <- subsetByOverlaps(CBpeaks, NBpeaks, minoverlap=1, ignore.strand=T) 
nbol = olRanges(NBpeaks, CBpeaks, ignore.strand=T)
cbol = olRanges(CBpeaks, NBpeaks, ignore.strand=T)

NB.unique = unique(NBpeaks[!(NBpeaks$peakName.NB) %in% (nbol$peakName.NB),])
NB.unique$peakName.CB <- NA
#NB.unique$foldChange.CB <- NA
#NB.unique$neg.logPv.CB <- NA
#NB.unique$relSummitPos.CB <- NA
mcols(NB.unique)$peakset = "NB.unique"

#cbol.any <- subsetByOverlaps(CBpeaks, NBpeaks, minoverlap=1, maxgap = 0, ignore.strand=T)

CB.unique <- unique(CBpeaks[!(CBpeaks$peakName.CB) %in% (cbol.any$peakName.CB),])
mcols(CB.unique)$peakset = "CB.unique"

CB.unique$peakName.NB <- NA
#CB.unique$foldChange.NB <- NA
#CB.unique$neg.logPv.NB <- NA
#CB.unique$relSummitPos.NB <- NA
mcols(CB.unique)$peakset = "CB.unique"

#for peaks with >75% overlap, calculate a new common peak as the intersection between the overlapping peaks

##75 percent overlaps, relative to nb, and relative to cb
##this is the same set, but the ranges are relative to NB or CB
cb75 =nbol[nbol$OLpercQ >= 75 | nbol$OLpercS >= 75]
nb75 = cbol[cbol$OLpercQ >= 75 | cbol$OLpercS >= 75]

seqnames(nb75) == seqnames(cb75)

mc.nb = mcols(nb75)
mc.cb = mcols(cb75)
mc.common <- cbind(mc.nb, mc.cb)
mc.common$peakName <- paste(mc.common$peakName.NB, mc.common$peakName.CB, sep="-")
names(mc.common)
mc.common <- mc.common[,c("peakName", "peakName.CB","peakName.NB", "IDRscore")]
mcols(nb75) <- mc.common
mcols(cb75) <- mc.common
names(nb75) <- mc.common$peakName
names(cb75) <- mc.common$peakName



all.75 <- (pintersect(nb75, cb75))
mcols(all.75) <- mc.common
all.75$peakset <- "CommonPeaks"

all.75 <- unique(all.75)



##Peaks that overlap leess than 75%
nbCut = GRangesList(nbol[nbol$OLpercQ < 75 & nbol$OLpercS < 75])
cbCut = GRangesList(cbol[cbol$OLpercQ < 75 & cbol$OLpercS < 75])



#for peaks with <75% overlap, extract protion with no overlap and create new unique peak

#nbcbCut = GRangesList(nbol[nbol$OLpercQ < 75])

#subtract overlap from CBpeaks for each peak with <75%
NBcut = unlist(psetdiff(nbCut, GRangesList(CBpeaks)))
CBcut <- unlist(psetdiff(cbCut, GRangesList(NBpeaks)))
#take these ranges from NBpeaks for metadata
NBcut.peaks <- subsetByOverlaps(NBpeaks, NBcut)
NBcut.peaks <- subsetByOverlaps(NBcut.peaks, NBcut)


median(width(NBcut.peaks))
NBcut.peaks$peakName.CB <- NA
#NBcut.peaks$foldChange.CB <- NA
#NBcut.peaks$neg.logPv.CB <- NA
#NBcut.peaks$relSummitPos.CB <- NA
mcols(NBcut.peaks)$peakset = "NB.cut"

##same for CB
CBcut = unlist(psetdiff(cbCut, GRangesList(NBpeaks)))





CBcut.peaks <- subsetByOverlaps(CBpeaks, CBcut)
CBcut.peaks <- subsetByOverlaps(CBcut.peaks, CBcut)
CBcut.peaks$peakName.NB <- NA
#CBcut.peaks$foldChange.NB <- NA
#CBcut.peaks$neg.logPv.NB <- NA
#CBcut.peaks$relSummitPos.NB <- NA
mcols(CBcut.peaks)$peakset = "CB.cut"






##get the intersection, for each range, for peaks with >75% overlap
##add this as commmon peak


save(Atlas.NBCB.common, file="atlas/Atlas.NBCB.common.Rdata")


##non-overlapping portions



#reduce(Atlas.NBCB.common)
#names(Atlas.NBCB.common)= paste(seqnames(Atlas.NBCB.common), start(Atlas.NBCB.common), "common", sep=":")

#Atlas.NBCB.common = reduce(Atlas.NBCB.common)

NB.unique$peakName <- NB.unique$peakName.NB
CB.unique$peakName <- CB.unique$peakName.CB
CBcut.peaks$peakName <-NA
NBcut.peaks$peakName <- NA   


names(mcols(NB.unique))
names(mcols(all.75))

Atlas.NBCB.idr =  unique(c(all.75[,c("peakset")], 
                       unique(NB.unique[,c("peakset")]), 
                       unique(CB.unique[,c("peakset")]),
                       unique(CBcut.peaks[,c("peakset")]),
                       unique(NBcut.peaks[,c("peakset")])))

Atlas.NBCB.idr[duplicated(Atlas.NBCB.idr),]

#####################################
##eliminate dups (overlaps existing)
##and maintain mcols
#######################
# m <- Atlas.NBCB.idr
# m2 <- reduce(m)
# length(m) == length(m2)
# IDX <- findOverlaps(m, m2)
# IDX2 <- IDX[which(!duplicated(subjectHits(IDX))),] #Just assign things once
# m3 <- m2[subjectHits(IDX2),]
# mcols(m3) <- mcols(m)[queryHits(IDX2),]
# 
# Atlas.NBCB.idr <- m3

save(Atlas.NBCB.idr, file="Atlas.NBCB.idr.Rdata")
save(Atlas.NBCB.idr, file="atlas/Atlas.NBCB.idr.Rdata")


#####ANNOTATE############################

##annotate atlas
library(ChIPpeakAnno)
library(ChIPseeker)
library(ChIPpeakAnno); library(GenomicFeatures)
data(TSS.human.NCBI36)
data("TSS.human.GRCh37")

#Atlas.NBCB.anno = endoapply(Atlas.NBCB, annotatePeakInBatch, AnnotationData=TSS.human.GRCh37)
require(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb = TxDb.Hsapiens.UCSC.hg19.knownGene
ge <- genes(txdb, columns=c("tx_name", "gene_id", "tx_type"))

Atlas.NBCB.idr.anno = annotatePeakInBatch(Atlas.NBCB.idr, FeatureLocForDistance = "TSS", select = "first",  AnnotationData=genes(txdb))
Atlas.NBCB.idr.anno <- addGeneIDs(Atlas.NBCB.idr.anno, silence=TRUE, orgAnn="org.Hs.eg.db", feature_id_type="entrez_id", IDs2Add=c("symbol"))

Atlas.NBCB.idr.anno$peakset <- Atlas.NBCB.idr$peakset
save(Atlas.NBCB.idr.anno, file="Atlas.NBCB.idr.anno.Rdata")
save(Atlas.NBCB.idr.anno, file="atlas/Atlas.NBCB.idr.anno.Rdata")

split(Atlas.NBCB.idr.anno, mcols(Atlas.NBCB.idr.anno)$peakset)

addGeneIDs(Atlas.NBCB.anno, IDs2Add = c("symbol"))




mcols(Atlas.NBCB.idr.anno)$
Atlas.NBCB.idr.anno[(Atlas.NBCB.idr.anno$distancetoFeature >3000 |Atlas.NBCB.idr.anno$distancetoFeature < -3000), ]


Atlas.NBCB.anno <- Atlas.NBCB.idr.anno

Atlas.NBCB.anno$peakset.coll <- Atlas.NBCB.anno$peakset

Atlas.NBCB.anno$peakset.coll[Atlas.NBCB.anno$peakset=="CB.cut"] <- "CB.unique"


Atlas.NBCB.anno$peakset.coll[Atlas.NBCB.anno$peakset=="NB.cut"] <- "NB.unique"

atlas.exp <- Atlas.NBCB.anno
Atlas.NBCB.idr.anno <- Atlas.NBCB.anno

names(atlas.exp) <- paste(names(atlas.exp), atlas.exp$peakset.coll, sep="-")

atlas.exp.deep <- atlas.exp[,"peakset.coll"]
#atlas.exp <- split(atlas.exp, mcols(atlas.exp)$peakset.coll)

export.bed(atlas.exp.deep, "Deeptools/atlas.export.bed")
export.bed(atlas.exp.deep, "AtlasFinal/atlas.export.bed")




####r
an= annotatePeakInBatch(Atlas.NBCB.idr, AnnotationData=genes(txdb))


##################################################################################
############################  PEAK COUNTS ####################################
##
samps = read.csv("samples.atac.csv", header=TRUE, stringsAsFactors = FALSE)
#samples_NBCB <- summarizeOverlaps(Atlas.idr, samps$bamReads, singleEnd=FALSE, ignore.strand=TRUE, fragments=FALSE, mode = "IntersectionNotEmpty")

Atlas.idr.l <- split(Atlas.idr, mcols(Atlas.idr)$peak)
samps = read.csv("samples.atac.T.csv", header=TRUE, stringsAsFactors = FALSE)

samples_NBCB <- samples_NBCB_atlas_list <- summarizeOverlaps(Atlas.idr.l, samps$bamReads, singleEnd=FALSE, ignore.strand=TRUE, fragments=FALSE, mode = "IntersectionNotEmpty")
save(samples_NBCB, file="samples_NBCB.Rdata")
# samples_NBCB_withfrags <- summarizeOverlaps(Atlas.idr.l, samps$bamReads, singleEnd=FALSE, ignore.strand=TRUE, fragments=TRUE, mode = "IntersectionNotEmpty")
# save(samples_NBCB_withfrags, file="samples_NBCB_withfrags.Rdata")
# samples_NBCB <-  samples_NBCB_OLatlas.list


#samples_NBCB_withfrags <- updateObject(samples_NBCB_withfrags)
sampsB <- samps[1:4,]
sampsB$Cell.type <- as.factor(sampsB$Cell.type)
sampsB$Cell.type <- relevel(sampsB$Cell.type, ref = "NB")
colData(samples_NBCB) <- DataFrame(sampsB)
#colData(samples_NBCB) <- colData(atacNB.exp.SE)
#colData(samples_NBCB_withfrags) <- colData(atacNB.exp.SE)

keep <- rowSums(assays(samples_NBCB_withfrags)$counts>10) >= 2
samples_NBCB_withfrags[keep,]

atacT100 <- atacT100[keep,]
samples_NBCB_withfrags

library(DESeq2)
ddsAtac <- DESeqDataSet(samples_NBCB, design = ~ Factor + Cell.type)
ddsAtac <- estimateSizeFactors(ddsAtac)
ddsAtac <- DESeq(ddsAtac)
results(ddsAtac)

atacNB.exp.SE <- SummarizedExperiment(assays = assays(ddsAtac)$counts, rowRanges=rowRanges(ddsAtac, drop=T), colData= colData(ddsAtac)) 
ddsAtac <- nor
#atacNB.exp.frags.SE <- SummarizedExperiment(assays = assays(ddsAtacf)$counts, rowRanges=rowRanges(ddsAtacf, drop=T), colData= colData(ddsAtacf)) 
###################UPDATE COUNTS ON ATLAS from atacNB.exp.SE summm exp with per-peaks counts that were summ overlaps on grangeslist#################

atlas <- unlist(rowRanges(atacNB.exp.SE))
Atlas.idr <- atlas 

#######################
######################


res <- resAtac <- results(ddsAtac)


plotMA(ddsAtac)
  
dds <- ddsAtac 


baseMeanNBrld = rowMeans(assay(rlog(dds, blind = FALSE))[,dds$BCell.type == "NB"])

baseMeanCBrld = rowMeans(assay(rlog(dds, blind = FALSE))[,ddsAtac$BCell.type == "CB"])

baseMeanNB = rowMeans(assay(dds, blind = FALSE)[,ddsAtac$BCell.type == "NB"])

baseMeanCB = rowMeans(assay(dds, blind = FALSE)[,ddsAtac$BCell.type == "CB"])



resAtac = cbind(baseMeanNB , baseMeanCB, baseMeanNBrld, baseMeanCBrld, as.data.frame(resAtac))
mcols(atacNB.exp.SE)$log2FoldChange <- resAtac$log2FoldChange
mcols(atacNB.exp.SE)$baseMeanNB <- resAtac$baseMeanNB
mcols(atacNB.exp.SE)$baseMeanCB <- resAtac$baseMeanCB
mcols(atacNB.exp.SE)$baseMeanNBrld <- resAtac$baseMeanNBrld
mcols(atacNB.exp.SE)$baseMeanCBrld <- resAtac$baseMeanCBrld




####

mcols(Atlas.idr)$symbol <- as.character(mcols(Atlas.idr)$symbol)




rownames(atacNB.exp.SE) <- names(Atlas.idr)

mc <- mcols(Atlas.idr)

ms <- mcols(atacNB.exp.SE)


mcols(atacNB.exp.SE) <- cbind(ms, mc)

save(atacNB.exp.SE, file="atacNB.exp.SE.Rdata")

#already done
mcols(Atlas.idr) <-  cbind(ms, mc)

#mcols(Atlas.idr) <- mcols(Atlas.idr)[,c(seq(37,66))]

save(Atlas.idr, file="Atlas.idr.Rdata")



##Done makginsummaried experiment on atlas ####
###############################################

#overlaps <- findOverlaps(Atlas.NBCB.idr.anno, atacNB.exp.SE, ignore.strand=T)


###export for deeptools


#### ATLAS PEAKS FOR CHIPSEQ REGIONS

#NB.CB.enhacers.unlist <- unlist(NB.CB.enhacers)

# overlaps <- findOverlaps(Atlas.idr, chipSeq_peaks.u, ignore.strand=T)
# 
# #create the match
# 
# 
# match_hit <- data.frame(names(Atlas.idr)[queryHits(overlaps)],
#                         names(chipSeq_peaks.u)[subjectHits(overlaps)],
#                         stringsAsFactors=F
# )
# 
# names(match_hit) <- c('query','subject')
# chp <-  mcols(chipSeq_peaks.u[match_hit$subject,])[,c("chpSeq")]
# atlas.chip<- Atlas.idr[match_hit$query,]
# mc.combined <- cbind(mcols(atlas.chip), data.frame(chp))
# atlas.chip.peaks <- granges(Atlas.NBCB.anno[match_hit$query,])
# mcols(atlas.chip.peaks) <- mc.combined
# split(atlas.chip.peaks, mcols(atlas.chip.peaks)$chp)
# 
# 
# 
# atlas.chip.peaks.list <- split(atlas.chip.peaks, mcols(atlas.chip.peaks)$chp)
# atlas.chip.peaks.list <- unique(atlas.chip.peaks.list)
# save(atlas.chip.peaks.list, file="atlas.chp.peaks.Rdata")
# #cal unique for each chhipseq exp
# atlas.chip.peaks <- unlist(atlas.chip.peaks.list)
# 
# 
# aebed <- atlas.chip.peaks[, c("chp")]
# 
# names(aebed) <- atlas.chip.peaks$chp
# aebed[aebed$chp== "CB_CTCF_T15F2_Peaks",]
# export.bed(aebed, "../MethSeq/atlas.chipseq.bed")
# 
# 
# ctcf <- aebed[aebed$chp== "CB_CTCF_T15F2_Peaks",]
# export.bed(ctcf, "../MethSeq/ctcf.bed")
# 
# bcor <- aebed[aebed$chp=="GCB_BCOR_HIseqt15f2c1",]
# export.bed(bcor, "../MethSeq/bcor.bed")
# 

##################ANNOTATION ##################
library(ChIPpeakAnno)
if(require(TxDb.Hsapiens.UCSC.hg19.knownGene)){
  aCR.allAtlas<-assignChromosomeRegion(Atlas.idr, nucleotideLevel=FALSE, proximal.promoter.cutoff=2000L,
                                       immediate.downstream.cutoff=2000L,
                              precedence=c("Promoters", "immediateDownstream",
                                           "fiveUTRs", "threeUTRs", 
                                           "Exons", "Introns"), 
                              TxDb=TxDb.Hsapiens.UCSC.hg19.knownGene)
  barplot(aCR.allAtlas$percentage)
}


##for just sign peaks

res[res$padj < 0.1 & (res$log2FoldChange >2 | res$log2FoldChange < -2),]
Atlas.idr.sig <- Atlas.idr[res$padj < 0.1 & (res$log2FoldChange >2 | res$log2FoldChange < -2),] 



if(require(TxDb.Hsapiens.UCSC.hg19.knownGene)){
    aCR.sigAtlas<-assignChromosomeRegion(Atlas.idr.sig, nucleotideLevel=FALSE, proximal.promoter.cutoff=2000L,
                                         immediate.downstream.cutoff=2000L,
                                         precedence=c("Promoters", "immediateDownstream", 
                                                      "fiveUTRs", "threeUTRs", 
                                                      "Exons", "Introns"), 
                                         TxDb=TxDb.Hsapiens.UCSC.hg19.knownGene)
    barplot(aCR.sigAtlas$percentage)
}

aCR.allAtlas <-  as.data.frame(aCR.allAtlas)
aCR.allAtlas$feature <- rownames(aCR.allAtlas)
aCR.allAtlas$peaks <- "allPeaks"

aCR.sigAtlas <- as.data.frame(aCR.sigAtlas)
aCR.allAtlas$feature <- rownames(aCR.allAtlas)
aCR.sigAtlas$feature <- rownames(aCR.sigAtlas)

aCR.sigAtlas$peaks <- "diffPeaks"


aCR <- rbind(as.data.frame(aCR.allAtlas), as.data.frame(aCR.sigAtlas))
#aCR <- data.frame(feature = rbind(rownames(aCR), rownames(aCR)), all = aCR.allAtlas$percentage, diffPeaks = aCR.sigAtlas$percentage)
#aCR$feature <-(rownames(aCR))
library(ggplot)

library(ggthemes)
ggplot(aCR, aes(peaks, percentage, fill=feature)) + geom_bar(stat = "identity") + scale_fill_hc() + theme_few()


aCR$percentage
save(aCR, file="aCR.chromRegions.Rdata")
source("https://bioconductor.org/biocLite.R")
library("ChIPseeker")


promoter <- getPromoters(TxDb=txdb, upstream=3000, downstream=3000)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
require(clusterProfiler)
plotAvgProf(tagMatrix, xlim=c(-3000, 3000), 
            xlab="Genomic Region (5'->3')", ylab = "Read Count Frequency")


atlas.ann <- annotatePeak(Atlas.idr, tssRegion=c(-3000, 3000), 
             TxDb=txdb, annoDb="org.Hs.eg.db")

granges(atlas.ann)
plotAnnoPie(atlas.ann)
plotAnnoBar(atlas.ann)
vennpie(atlas.ann)
plotDistToTSS(atlas.ann)
#could proivde bedgraphs
plotAvgProf2(c("NB.idr.bed", "CB.idr.bed"))

library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb = TxDb.Hsapiens.UCSC.hg19.knownGene
anno <- annotatePeak(Atlas.NBCB.idr.anno, tssRegion = c(-2000, 2000), TxDb = txdb, annoDb="org.Hs.eg.db")



mcols(Atlas.idr)$peakset <- mcols(Atlas.NBCB.idr.anno)$peakset
mcols(Atlas.idr)$peak <- mcols(Atlas.NBCB.idr.anno)$peak
save(Atlas.idr, file="Atlas.idr")



#####################Expression Data ###########################
##expdata in Atlas
load("/Volumes/AD/Projects/ATAC/Normal_B_Cells/ddsNB.rna-seq-directional.Rdata")
resNBCexp <- results(ddsNB)

baseMeanNB.exp = rowMeans(assay(rlog(ddsNB, blind = FALSE))[,ddsNB$BCell.type == "NB"])

baseMeanCB.exp = rowMeans(assay(rlog(ddsNB, blind = FALSE))[,ddsNB$BCell.type == "CB"])

resNBCexp = results(ddsNB)

resNBCexp= cbind(baseMeanNB.exp , baseMeanCB.exp, as.data.frame(resNBCexp))

resNBCexp$symbol <- rownames(resNBCexp)

resNBCexp$log2FC.exp <- resNBCexp$log2FoldChange
resNBCexp[is.na(resNBCexp$log2FC.exp ),"log2FC.exp"] <- 0
resNBCexp[is.na(resNBCexp$padj),"padj.TF"] <- 1


nb.match <- resNBCexp %>% select(symbol, log2FC.exp)
nb.match <- na.omit(nb.match)


mc <- as.data.frame(mcols(Atlas.idr))

mc$symbol <- as.character(mc$symbol)
mcn <- left_join(mc, nb.match, by="symbol")

length(rownames(mcn))
length(rownames(mc))

mcols(Atlas.idr) <- mcn
##############################
#reintegrate with chip-HiC data

##can use the original summarized experiment for atlas olap, since this is baed on reads of atlas

###do findoverlaps on atlas
# overlaps <- findOverlaps(Atlas.idr, chp, ignore.strand=T)
# 
# match_hit <- data.frame(names(Atlas.idr)[queryHits(overlaps)],
#                         names(chp)[subjectHits(overlaps)],
#                         stringsAsFactors=F
# )
# 
# names(match_hit) <- c('query','subject')
# n[seq(38,42)]
# n <- colnames(mcols(Atlas.idr))
# nn <- n[!n %in% colnames(mcols(chp))]
# newmcM <-  mcols(Atlas.idr[match_hit$query,])[,c("peak", "log2FC.exp", "value")]
# chp.ov <- chp[match_hit$subject,]
# 
# mc.combined <- cbind(mcols(chp.ov), newmcM)
# 
# mcols(chp.ov) <- mc.combined 
# 
# chip.HiC.atlas.peaks <- split(chp.ov, mcols(chp.ov)$chpSeq, drop = F)
# 
# save(chip.HiC.atlas.peaks, file="chip.HiC.atlas.peaks.Rdata")
# save(chip.HiC.atlas.peaks, file="AtlasFinal/chip.HiC.atlas.peaks.Rdata")
# 
# save(Atlas.idr, file="AtlasFinal/Atlas.ird.Rdata")

###DONE###

# overlaps <- findOverlaps(atacNB.exp.SE, chp, ignore.strand=T)
# 
# match_hit <- data.frame(rownames(atacNB.exp.SE)[queryHits(overlaps)],
#                         names(chp)[subjectHits(overlaps)],
#                         stringsAsFactors=F
# )
# 
# names(match_hit) <- c('query','subject')
# n[seq(38,42)]
# n <- colnames(mcols(atacNB.exp.SE))
# newmcM <-  mcols(atacNB.exp.SE[match_hit$query,])[,c(n[seq(38,42)], "peakset.coll", "value")]
# chp.ov <- chp[match_hit$subject,]
# 
# mc.combined <- cbind(mcols(chp.ov), newmcM)
# 
# mcols(chp.ov) <- mc.combined 
# 
# chip.HiC.atlas.peaks <- split(chp.ov, mcols(chp.ov)$chpSeq, drop = F)
# 
# save(chip.HiC.atlas.peaks, file="chip.HiC.atlas.peaks.Rdata")
# ###DONE###

###just add expression...



Atlas.idr$peakset[Atlas.idr$peakset=="CB.cut"] <- "CB.unique"

#####export
Atlas.idr.exp <- Atlas.idr
names(Atlas.idr.exp) <- mcols(Atlas.idr.exp)$peakset
save(Atlas.idr, file="AtlasFinal/Atlas.idr.Rdata")
names(Atlas.idr.exp) <- mcols(Atlas.idr.exp)$peakset
save(Atlas.idr, file="Atlas.idr.Rdata")
export.bed(Atlas.idr.exp, "Atlas.idr.bed")

export.bed(atlas.idr.exp, "atlas.exp.bed")
#export with +/- 2.5kb


width(Atlas.idr)
Atlas.idr.5kb <- resize(Atlas.idr, 5000, fix = "center")
Atlas.idr.5kb <- trim(Atlas.idr.5kb)
export.bed(Atlas.idr.5kb, "Atlas.idr.5kb.bed")


#########################Atlas x chipseq#################################
#relative to atlas, output is grangeslist, one granges for each chip-seq
atlas <- na.omit(Atlas.idr)
names(atlas) <- NULL

chipSeq_peaks

mcols(atlas)[is.na(mcols(atlas)),]
Atlas.idr.chip <- GRangesList(mclapply(chipSeq_peaks, function(x) {
  subsetByOverlaps(Atlas.idr, x, ignore.strand=T)
},  mc.cores = getOption("mc.cores", 2L)))
Atlas.idr.chip <- unique(Atlas.idr.chip)
#save(Atlas.idr.chip, file="AtlasFinal/Atlas.idr.chip.Rdata")
save(Atlas.idr.chip, file="Atlas.idr.chip.Rdata")

chip.df <- tbl_df(as.data.frame(Atlas.idr.chip, group_name.as.factor=TRUE))


#############################

x <- lapply(Atlas.idr.chip, names)
Atlas.idr.chip.u <- unlist(Atlas.idr.chip, use.names = T)
relist(Atlas.idr.chip.u)
mcols(Atlas.idr.chip.u)$TF <- names(Atlas.idr.chip.u)
save(Atlas.idr.chip.u, file="AtlasFinal/Atlas.idr.chip.u.Rdata")








##clean up mcols and resave
# mcols(Atlas.idr) <- (mcols(Atlas.idr))[c(seq(37,65),68)]
# 
# 
# mcols(Atlas.idr)$log2FC.exp = mcols(Atlas.idr)$log2FC.exp.y
# 
# 



##to plotting
res <- data.frame(mcols(Atlas.idr.chip.u))

resC <- res

resC$TF <- as.factor(resC$TF)
resC$TF <- reorder(resC$TF, resC$log2FoldChange, FUN=median)

##Got to chipseqAtlas.R file, line ~240 for plots


##Atlas x enhancers
atlas <- Atlas.idr
Atlas.idr.enhancers <- GRangesList(mclapply(NB.CB.enhacers, function(x) {
  subsetByOverlaps(atlas, x, ignore.strand=T)
},  mc.cores = getOption("mc.cores", 2L)))



Atlas.idr.enhancers <- GRangesList(mclapply(NB.CB.enhacers, function(x) {
  subsetByOverlaps(Atlas.idr, x, ignore.strand=T)
},  mc.cores = getOption("mc.cores", 2L)))



en.df <- as.data.frame(Atlas.idr.enhancers, group_name.as.factor = TRUE)


en.df$var <- rowVars(log2(en.df[,38:41]+1))

var.region <- ddply(en.df, .(group_name), summarize,
                     var = unique(var),
                     ecdf = ecdf(var)(unique(var)))     

mat <- en.df[,38:41]
qmat <- normalize.quantiles.robust(as.matrix(mat))
hist(en.df[,38:41])
hist(qmat)
en.df[,38:41]<- qmat 

ggplot(var.region, aes(var, ecdf, color = group_name)) + geom_step() + scale_color_brewer(palette = "Set1") + theme_few()



###variance distribution on chipseq regions
chip.df <- tbl_df(as.data.frame(Atlas.idr.chip, group_name.as.factor=TRUE))
chip.df$chip.grp <- NA

grpCB <- c("CB_CTCF_T15F2_Peaks", "B_IKZF1_T15F2_Peaks", "CB051710_PU1_T15F2_Peaks", "CB22_HSF1_T15F2_targets", "GCB_BCL6peaks_t15f2f158c1", "GCB_BCOR_HIseqt15f2c1", "GCB_SMRT_HIseqt15f2c1", 
           "GCB_SOX9_hg18_T10F2_BR1overlapBR2.bed", "GCB_SMRT_HIseqt15f2c1", "GCB_SOX9_hg18_T10F2_BR1overlapBR2.bed")

grpLY1 <- c("Ly1_BCL6peaks_r1r3_hg19.bed", "Ly1_BCORab185_t15f2f158c1", "LY1_CTCF_T15F2_Peaks", "Ly1_NCORpeaks_t15f2f158c1", "LY1_P300_KC_T15F2_Peaks", "Ly1_SMRTpeaks_t15f2c1", "LY1100518_PAX5_T15F2_Peaks")
grpLY7 <- c("LY7_BCL6rep2_t15f2f158c1", "LY7_BCORpeaks_t15f2f158c1", "LY7_SMRTpeaks_t15f2f158c1")
grpNB <- c( "NB_CTCF_T15F2_Peaks", "NB051710_PU1_T15F2_Peaks", "NB20_22_HSF1_T15F2_targets")

chip.df$chip.grp[chip.df$group_name %in% grpLY1] <- "LY1"
chip.df$chip.grp[chip.df$group_name %in% grpCB] <- "CB"
chip.df$chip.grp[chip.df$group_name %in% grpNB] <- "NB"


chip.df$var <- rowVars(log2(chip.df[,43:46]+1))
chip.df$var <- rowVars(log2(chip.df[,c("NB.58.qcnt", "CB.58.qcnt", "NB.59.qcnt", "CB.59.qcnt")]))
library(plyr)
LY1.var.chip <- ddply(chip.df[chip.df$chip.grp %in% "LY1",], .(group_name), summarize,
                    var = unique(var),
                    ecdf = (1- ecdf(var)(unique(var))))

ggplot(LY1.var.chip, aes(var, (ecdf), color = group_name)) + geom_step() + scale_color_brewer(palette = "Set1") + theme_few() + xlim(c(0, 5))


CB.var.chip <- ddply(chip.df[chip.df$chip.grp %in% "CB",], .(group_name), summarize,
                      var = unique(var),
                      ecdf = (1- ecdf(var)(unique(var))))

ggplot(NB.var.chip, aes(var, (ecdf), color = group_name)) + geom_step() + scale_color_brewer(palette = "Set1") + theme_few() + xlim(c(0, 5))


NB.var.chip <- ddply(chip.df[chip.df$chip.grp %in% "NB",], .(group_name), summarize,
                     var = unique(var),
                     ecdf = (1- ecdf(var)(unique(var))))

ggplot(NB.var.chip, aes(var, (ecdf), color = group_name)) + geom_step() + scale_color_brewer(palette = "Set1") + theme_few() + xlim(c(0, 5))

#####Annotoation figures for peak atlas


library(ChIPseeker)







##############MOTIFS######################
load("TF.ranges.u.Rdata") ##motif ranges from homer
TF.ranges.u

atlas.motifs <- GRangesList(mclapply(TF.ranges.u, function(x) {
    subsetByOverlaps(atlas, x, ignore.strand=T)
},  mc.cores = getOption("mc.cores", 2L)))

save(atlas.motifs, file="atlas.motifs.Rdata")


        