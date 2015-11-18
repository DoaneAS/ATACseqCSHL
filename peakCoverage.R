library(GenomicRanges)
library(systemPipeR)
library(GenomicAlignments)
args <- systemArgs(sysma="param/count_rangesets.param", mytargets="targets_macs.txt")
args_bam <- systemArgs(sysma=NULL, mytargets="targets_bam.txt")
bfl <- BamFileList(outpaths(args_bam), yieldSize=50000, index=character())
countDFnames <- countRangeset(bfl, args, mode="Union", ignore.strand=TRUE)
writeTargetsout(x=args, file="targets_countDF.txt", overwrite=TRUE)

library(rtracklayer)
samps = read.csv("samples.atac.T.csv", header=TRUE, stringsAsFactors = FALSE)

bfl = samps$bamReads[1]

export.bed(reduce(Atlas.NBCB.union.anno[1:100,]), "atlas/test.bed")

peaks <- import.bed("Atlas.idr.5kb.bed")
names(peaks) <- paste0(as.character(seqnames(peaks)), "_", start(peaks), "-", end(peaks))
peaks <- split(peaks, names(peaks))

countDF.s1 <- summarizeOverlaps(peaks, bfl)
countDF <- assays(countDF)$counts
countDF
write.table(countDF, "sample_2.atlas.cov.Rdata", col.names=NA, quote=FALSE, sep="\t")


bfl = samps$bamReads[3]
countDF <- summarizeOverlaps(peaks, bfl)
countDF <- assays(countDF)$counts
write.table(countDF, "sample_3.atlas.cov.Rdata", col.names=NA, quote=FALSE, sep="\t")

bfl = samps$bamReads[4]
countDF <- summarizeOverlaps(peaks, bfl)
countDF <- assays(countDF)$counts
write.table(countDF, "sample_4.atlas.cov.Rdata", col.names=NA, quote=FALSE, sep="\t")


colData(atacNB.exp.SE) <- colData(dds)

colData(atacNB.exp.SE)$BCell.typ <- colData(dds)$Cell.type
colData(atacNB.exp.SE)$PatientID <- colData(atacNB.exp.SE)$Factor
# sample_1.atlas.counts <- summarizeOverlaps(peaks, samps$bamReads[1])
# sample_2.atlas.counts <- summarizeOverlaps(peaks, samps$bamReads[2])
# sample_3.atlas.counts <- summarizeOverlaps(peaks, samps$bamReads[3])
# sample_4.atlas.counts <- summarizeOverlaps(peaks, samps$bamReads[4])
# 
# newCountDataSet(assay(sample_1.atlas.counts), rownames(colData(sample_1.atlas.counts)))

####COUNTS DATA ON PEAKS#####

##count using summarize overlaps###

length(unique(Atlas.idr$peak))
atlas.list <- split(Atlas.idr, mcols(Atlas.idr)$peak)
samples_NBCB_OLatlas.list <- summarizeOverlaps(atlas.list, samps$bamReads, singleEnd=FALSE, ignore.strand=TRUE, fragments=FALSE)

save(samples_NBCB, file="samples_NBCB.sumEx.atac.idr.Rdata")
load("samples_NBCB.sumEx.atac.Rdata")
NBCBpheno = data.frame(SampleID = samps$SampleID, BCell.type=samps$Condition, PatientID = samps$Factor, bamReads=samps$bamReads)

NBCBpheno$BCell.type <- relevel(NBCBpheno$BCell.type, "NB")
rownames(NBCBpheno) <- NBCBpheno$SampleID
colData(samples_NBCB)= colData(atacNB.exp.SE) #DataFrame(NBCBpheno)

colData( samples_NBCB_OLatlas.list) <- colData(atacNB.exp.SE)
#############

#for DESeq
#NBCB_counts <- newCountDataSet(assay(samples_NBCB), rownames(colData(samples_NBCB)))

head(colData(samples_NBCB))
################################
###########for DESEq2       ####
#################################
library(DESeq2)
NBCB_pcounts <- DESeqDataSet(samples_NBCB_OLatlas.list, design = ~ PatientID + BCell.type)
NBCB_pcounts
ddsNBCBatac <- estimateSizeFactors(NBCB_pcounts)
ddsNBCBatac = DESeq(ddsNBCBatac)
ddsNBCBatac
res.atac <- results(ddsNBCBatac)
res.atac

baseMeanNBatac = rowMeans(assay(rlog(ddsNBCBatac))[,ddsNBCBatac$BCell.type == "NB"])
baseMeanCBatac = rowMeans(assay(rlog(ddsNBCBatac))[,ddsNBCBatac$BCell.type == "CB"])
head(baseMeanCB)

res.atac = cbind(baseMeanNBatac, baseMeanCBatac, as.data.frame(res.atac))
res.atac = cbind(peakNumber=rownames(res.atac), as.data.frame(res.atac))

res.atac$covName = rownames(res.atac)

res.atac = data.frame(covName=res.atac$covName, baseMeanNBatac=res.atac$baseMeanNBatac, baseMeanCBatac=res.atac$baseMeanCBatac,log2FC_CBvNB_atac.counts=res.atac$log2FoldChange, padj.atac.counts=res.atac$padj, peakName=as.character(res.atac$symbol))

rownames(res.atac) = as.character(res.atac$covName)
#link to atlas
mcols(Atlas.NBCB.union.anno.sym)$covName = paste0(as.character(seqnames(Atlas.NBCB.union.anno.sym)), "_", start(Atlas.NBCB.union.anno.sym), "-", end(Atlas.NBCB.union.anno.sym))

df = as.data.frame(Atlas.NBCB.union.anno.sym)
df$covName <- as.character(df$covName)
mcols(Atlas.NBCB.union.anno.sym)$covName %in% rownames(res.atac)
  
library(dplyr)

NBCB.pcounts = le_join(df, res.atac, by="covName")
NBCB.pcounts


NBCB.exp.atac$log2FC_NBvCB_rna


####coverage on specific regions
##bivalent


summarizeOverlaps(peaks, samps$bamReads, singleEnd=FALSE, ignore.strand=TRUE, fragments=FALSE)

######################################
#BIVALENT#### coverage over all bivalent regions
##################################
library(rtracklayer)
library(GenomicAlignments)

GCBbiv = import.bed("/Volumes/AD/Projects/ATAC/MethSeq/bivalent/GCB_denovoBivalentPromoters.hg19.bed")
GCBbiv
GCBcont = import.bed("/Volumes/AD/Projects/ATAC/MethSeq/bivalent/GCB_K4mono_controlPromoters.hg19.bed")
mcols(GCBcont)$BivCont <- "Control"
mcols(GCBbiv)$BivCont <- "Biv"
GCBbivAcont = c(GCBbiv, GCBcont)

mcols(GCBbivAcont)$symbol = as.character(mcols(GCBbivAcont)$name)

names(GCBbivAcont) <- paste0(as.character(seqnames(GCBbivAcont)), "_", start(GCBbivAcont), "-", end(GCBbivAcont))
GCBbivAcont <- split(GCBbivAcont, names(GCBbivAcont))
GCBbivAcont[mcols(GCBbivAcont)$BivCont=="Biv",]

GCBbivAcont[which(mcols(GCBbivAcont)$BivCont=="Biv"),]

NBCB.atac.biv.cont = summarizeOverlaps(GCBbivAcont, samps$bamReads, singleEnd=FALSE, ignore.strand=TRUE, fragments=FALSE)



NBCBpheno = data.frame(SampleID = samps$SampleID, BCell.type=samps$Condition, PatientID = samps$Factor, bamReads=samps$bamReads)
NBCB.atac.biv.cont
NBCBpheno$BCell.type <- relevel(NBCBpheno$BCell.type, "NB")
rownames(NBCBpheno) <- NBCBpheno$SampleID

colData(NBCB.atac.biv.cont)= DataFrame(NBCBpheno)

rowRanges(NBCB.atac.biv.cont)

ddsBiv = DESeqDataSet(NBCB.atac.biv.cont, design = ~ PatientID + BCell.type)
ddsBiv = DESeq(ddsBiv)
rowRanges(ddsBiv)
rlog(ddsBiv)
resBiv = results(ddsBiv)
resBiv
ddsBiv
plotMA(resBiv)

resGr = results(ddsBiv, lfcThreshold=1, format="GRanges")
resGr
mcols(resGr)$symbol <- as.character(mcols(GCBbivAcont)$name)

mcols(resGr)$BivCont <- as.character(mcols(GCBbivAcont)$BivCont)

resGr.biv = resGr[resGr$BivCont=="Biv",]

resGr = as.data.frame(resGr)

head(resGr)
library(ggplot2)
ggplot(resGr, aes(x=BivCont, y=log2FoldChange)) + geom_point(position=position_jitter(width=.3, height=.48), alpha=0.3, size=1.5) + geom_boxplot(width=0.7, alpha=0, outlier.size=0, color="black") +
  theme(legend.position="none") +  scale_colour_brewer(palette="Set1") + theme_minimal() +
  ylab("CB vs NB\n log2FC") + theme(legend.position="none")

rld = rlog(ddsBiv)


rld.biv = rld[which(mcols(GCBbivAcont)$BivCont=="Biv"),]
rld.cont = rld[which(mcols(GCBbivAcont)$BivCont!="Biv"),]

(rld.biv) <- mcols(GCBbivAcont)$name[which(mcols(GCBbivAcont)$BivCont=="Biv")]

rld.df = (as.data.frame(assay(rld)))

rld.df$BivCont = as.factor(mcols(GCBbivAcont)$BivCont)
rld.df$BivCont = relevel(rld.df$BivCont, ref="Control")

rld.df[rld.df$BivCont=="Biv",]
results(ddsBiv, lfcThreshold=1)$log2FoldChange
rld.df$log2FC = results(ddsBiv, lfcThreshold=1)$log2FoldChange
rld.df$p.adj = results(ddsBiv, lfcThreshold=1)$p.adj



library(edgeR)


#########################################################
##bivalent- overlap and civerage on overlapping regions###
#########################################################

##make granges of updated 

Atlas.new <- Atlas.NBCB.union.anno.sym
names(Atlas.new) <- paste0(as.character(seqnames(Atlas.NBCB.union.anno.sym)), "_", start(Atlas.NBCB.union.anno.sym), "-", end(Atlas.NBCB.union.anno.sym))

mcols(Atlas.new) <- cbind( mcols$Atlas.new)




NBCB.exp.atac
  NBCB.exp.atac

source("https://bioconductor.org/biocLite.R")
biocLite("edgeR")



p <- colData(atacNB.exp.SE)
samples_NBCB <- updateObject(samples_NBCB)
cd <- data.frame(SampleID = c("NB_58", "CB-58", "NB-59", "CB-59"), BCell.type = p$BCell.type, PatientID = p$PatientID)
colData(samples_NBCB) <- DataFrame(cd)
colnames(samples_NBCB) <- c("NB_58", "CB-58", "NB-59", "CB-59")

save(samples_NBCB, file="samples_NBCB.Rdata")
dds <- DESeqDataSet(samples_NBCB, design = ~ PatientID + BCell.type)
dds$BCell.type <- relevel(dds$BCell.type, ref = "NB")
dds <- DESeq(dds)
res <- results(dds)


#######################
####################UPDATE ATAC SE ######################
library(DESeq2)
library(preprocessCore)
atacNB.exp.SE <- updateObject(atacNB.exp.SE)

save(atacNB.exp.SE, file="atacNB.exp.SE.Rdata")
mcols(atacNB.exp.SE)


mcols(Atlas.idr)$annotation
colData(atacNB.exp.SE)
names(assays(atacNB.exp.SE)) <- "counts"

dds <- DESeqDataSet(samples_NBCB, design = ~ PatientID + Cell.type)
dds <- estimateSizeFactors(dds)
dds$BCell.type <- relevel(dds$BCell.type, ref = "NB")
dds <- DESeq(dds)
res <- results(dds)
cnt <- counts(dds)
rld <-assay(rlog(dds, blind = FALSE))
q.cnt <- normalize.quantiles(cnt)

colnames(cnt) <- c("NB_58.cnt", "CB-58.cnt", "NB-59.cnt", "CB-59.cnt")
colnames(rld) <- c("NB-58.rld", "CB-58.rld", "NB-59.rld", "CB-59.rld")
colnames(q.cnt)  <- c("NB-58.qcnt", "CB-58.qcnt", "NB-59.qcnt", "CB-59.qcnt")


mcols(Atlas.idr)$padj <- res$padj
save(Atlas.idr, file="Atlas.idr.Rdata")

mc <- mcols(Atlas.idr)[,1:30]

colnames(mc)
mcols(Atlas.idr) <- cbind(mc, as.data.frame(cnt), as.data.frame(rld), as.data.frame(q.cnt))
atacNB.exp.SE <- SummarizedExperiment(assays = assays(dds)$counts, rowRanges=rowRanges(dds, drop=T), colData= colData(dds)) 
peaks <- tbl_df(as.data.frame( mcols(Atlas.idr)[,1:30]))


peaks$annotation <- as.character(peaks$annotation)

peaks$annotation

peaks$genomicR <- NA

peaks$genomicR[peaks$annotation %in% "Promoter (<=1kb)"] <- "Promoter1"


peaks$genomicR[peaks$annotation %in% "Promoter (1-2kb)"] <- "Promoter2"


peaks$genomicR[peaks$annotation %in% "Distal Intergenic"] <- "DistalIntergenic"

peaks <- cbind(peaks, as.data.frame(cnt), as.data.frame(q.cnt), as.data.frame(rld))
mcols(Atlas.idr) <-peaks
peaksR <- peaks %>% group_by(genomicR)
#add back to atlas

#mc <- mcols(Atlas.idr)
#qd <- as.data.frame(qmat)
#colnames(qd) <- c("NB1", "CB1", "NB2", "CB2")
#mcols(Atlas.idr) <- cbind(mcols(Atlas.idr), qd)



#mcols(Atlas.idr)[,40:43] <- rld


d <- plotCounts(dds, gene=which.min(res$padj), intgroup="BCell.type")#,returnData=TRUE)
library("ggplot2")
ggplot(d, aes(x=Bell.type, y=count)) +
  geom_point(position=position_jitter(w=0.1,h=0)) +
  scale_y_log10(breaks=c(25,100,400))
mat <- (peaksR[peaksR$genomicR %in% "Promnoter",32:35])

var(mat)

library(genefilter)
library(preprocessCore)



qmat <- normalize.quantiles(as.matrix(peaks[,36:39]))
peaks <- cbind(peaks, as.data.frame(qmat))
peaks$var <- rowVars(peaks[,32:35])
peaks$cnt.var <- rowVars(log2(peaks[,36:39]+1))
peaks$qvar <- rowVars(log2(peaks[,40:43]+1))
peaksR <- peaks %>% group_by(genomicR)
PromVar <- rowVars(peaksR[peaksR$genomicR %in% "Promoter",32:35])

plot(PromVar)

library(plyr)
cumfreqAll <- summarize(peaks, var = unique(var), 
                      ecdf = ecdf(var)(unique(var)) * length(var))
var.region <- ddply(peaks, .(genomicR), summarize,
                      cnt.var = unique(cnt.var),
                      ecdf = ecdf(cnt.var)(unique(cnt.var)))

qvar.region <- ddply(peaks, .(genomicR), summarize,
                    qvar = unique(qvar),
                    ecdf = ecdf(qvar)(unique(qvar)))

ggplot(var.region, aes(cnt.var, ecdf, color = genomicR)) + geom_step()

ggplot(qvar.region, aes(qvar, ecdf, color = genomicR)) + geom_step() + scale_color_brewer(palette = "Set1") + theme_few()


ggplot(df, aes(x, colour = g)) + stat_ecdf()

library(ggplot2)
library(ggthemes)


qplot(unique(peaks$var), ecdf(peaks$var)(unique(peaks$var)), geom='step')




###variance in enahncers regions


subsetByOverlaps(atacNB.exp.SE
                 

peaks




#############ATAC LOAD ####

library("rtracklayer")
atac.load <- import.bed("atlas.idr.12.5kb.bed")
atac.load
library(GenomicAlignments)

atac.load <- summarizeOverlaps(atac.load, samps$bamReads, singleEnd=FALSE, ignore.strand=TRUE, fragments=FALSE)


colData(atac.load) <- colData(atacNB.exp.SE)
atac.load.SE <- atac.load
save(atac.load.SE, file="atac.load.SE.Rdata")
cts <- assay(atac.load[,c(2)])
i <- (sort(cts, index.return=TRUE))
atac.load<- atac.load[i$ix,]
library(ChIPpeakAnno)
require(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb = TxDb.Hsapiens.UCSC.hg19.knownGene
ge <- genes(txdb, columns=c("tx_name", "gene_id", "tx_type"))


atac.load.anno = annotatePeakInBatch(rowRanges(atac.load), FeatureLocForDistance = "TSS", select = "first",  AnnotationData=genes(txdb))
atac.load.anno  <- addGeneIDs(atac.load.anno, silence=TRUE, orgAnn="org.Hs.eg.db", feature_id_type="entrez_id", IDs2Add=c("symbol"))
rowRanges(atac.load) <- atac.load.anno
library(DESeq2)
ddsAtLoad <- DESeqDataSet(atac.load, design = ~ PatientID + BCell.type) 
ddsAtLoad  <- estimateSizeFactors(ddsAtLoad )
ddsAtLoad <- DESeq(ddsAtLoad)

norm.quant.atac <- function(dds) {
    require(preprocessCore)
    require(GenomicRanges)
    if (class(dds)[1] == "DESeqDataSet") {
        nq <-normalize.quantiles(x=counts(dds, normalize=TRUE))
        assays(dds)$cts.qn <- log(nq[,]+1)/log(2)
        return(dds)
    } 
    else if (class(dds)[1] == "RangedSummarizedExperiment") {
        sumexp <- dds
        nq <-normalize.quantiles(x=assay(sumExp))
        assays(sumExp)$cts.qn <- log(nq[,]+1)/log(2)
        return(sumexp)
        
    }
}
ddsAtLoad <- norm.quant.atac(ddsAtLoad)
assays(atac.load)$cts.qn <- assays(ddsAtLoad)$cts.qn
res <- results(ddsAtLoad)
mcols(atac.load)$log2FoldChange <- res$log2FoldChange
mcols(atac.load)$padj <- res$padj




assays(atac.load)$cts.norm <- counts(ddsAtLoad, normalized=TRUE)
assays(atac.load)$rld <- assay(rlog(ddsAtLoad, blind=FALSE))

save(atac.load, file="atac.load.Rdata")

#mcols(atac.load)$log2.cts <- log2(assay(atac.load))

#p300.gr.anno.enh <- p300.gr.anno[abs(p300.gr.anno$distancetoFeature) > 2000,]



norm.quantiles.atac.old <- function(sumExp) {
    require(preprocessCore)
    require(GenomicRanges)
    nq <-normalize.quantiles(x=assay(sumExp))
    assays(sumExp)$cts.qn <- log(nq[,]+1)/log(2)
    return(sumExp) 
    
}



atac.load
# CB.mean <- rowMeans(assay(atac.load[,c(2,4)]))
# NB.mean <- rowMeans(assay(atac.load[,c(1,3)]))
# 
 plot(range01(assays(atac.load[,c(2)])$cts.qn))
 
plot(range01(assay(atac.load[,c(2)])))
mcols(atac.load)$CB.rank <- range01(assay(atac.load[,c(2)]))

plot(df$i
 
 
 hist(range01(assay(atac.load[,c(2)])))
 
 df <- as.data.frame(mcols(atac.load))
 
 df$i <- range01(assay(atac.load[,c(2)]))
df$CB.cts <- assay(atac.load[,c(2)])
x <- quantile(assay(atac.load[,c(2)]), probs=c(0.95))
df[df$CB.cts >x,]
highest.CB <- df[df$CB.cts >x,]


export.bed(rowRanges(atac.load[31000:32141,]), "atac.load.top.bed")

atac.top.load <- rowRanges(atac.load[31000:32141,])

findOverlaps(atac.top.load, p300SE)
subsetByOverlaps(atac.top.load, p300.SuperE)

#highest.CB <- highest.CB[(highest.CB$log2FoldChange>1 | highest.CB$log2FoldChange< (-1)),]
write.table(highest.CB[,c("log2FoldChange", "symbol", "CB.cts", "padj"),], "highestCB.txt")

plot(df$i)
text(df$i[32100:32141], labels=df$symbol[32100:32141], cex= 0.5, pos=3) 
    df
library(dplyr)
 df <- unique.data.frame(df) 
df <- tbl_df(df)
df
highest.FC <- df[(df$log2FoldChange>1 | df$log2FoldChange< (-1)),]



i <- sort(highest.FC$log2FoldChange, index.return=T)


write.table(highest.FC[i$ix,c("log2FoldChange", "symbol")], "highFG.CBload.txt")

atac.load.p300 <- subsetByOverlaps(atac.load, p300.enhancers.u)

plot(range01(assays(atac.load.p300[,c(2)])$cts.qn))

plot(range01(normalize.quantiles(assay(atac.load.p300[,c(4)]))))
plot(range01(assay(atac.load.p300[,c(4)])))

df <- as.data.frame(rowRanges(atac.load.p300))

df$i <- range01(assay(atac.load.p300[,c(4)]))
df$CB.cts <- assay(atac.load.p300[,c(4)])
x <- quantile(assay(atac.load.p300[,c(4)]), probs=c(0.9))
df[df$CB.cts >x,]
highest.CB <- df[df$CB.cts >x,]

highest.CB <- highest.CB[(highest.CB$log2FoldChange>1 | highest.CB$log2FoldChange< (-1)),]
write.table(highest.CB[,c("log2FoldChange", "symbol", "CB.cts"),], "highestCB.txt")




    )
# normalize.quantiles
# rowMeans(
# normTransform(
# range01 <- function(x){(x-min(x))/(max(x)-min(x))}
# 
# plot(range01(assays(atac.load)$cts.qn))
# plot(assays(atac.load)$cts.qn)
# 
# cb <- as.data.frame(range01(rowMeans(assays(atac.load[,c(2,4)])$cts.qn)))
# colnames(cb) <- "CB.score"
# 
# nb <- as.data.frame(range01(rowMeans(assays(atac.load[,c(1,3)])$cts.qn)))
# head(nb)
# colnames(nb) <- "NB.score"
# df.score <- cbind(nb, cb)
# df.score <- na.omit(df.score)
# rownames(df)
# df.score$i <- 1:nrow(df.score)
# df.score$i <- range01(df.score$i)  
# 
# head(df.score)
# df.score$NB.qn <- assays(atac.load[,c(1,3)])$cts.qn
# df.score$CB.qn <- assays(atac.load[,c(2,4)])$cts.qn
# df.score$NB.cts <- assays(atac.load[,c(1,3)])$counts
# df.score$CB.cts <- assays(atac.load[,c(2,4)])$counts
# length(df.score$CB.qn)
# plot(df.score$CB.qn, df.score$CB.score)



plot(df.score)
quantile(dx$score, probs=c(0.95))
quantile(mcols(atac.load)$log2cts, probs=c(0.90))
id <- mcols(p300.gr.anno.enh)$log2cts > 8.78
idd <- mcols(p300.gr.anno.enh)$log2cts < 8.78


p300.SuperE <- p300.gr.anno.enh[as.vector(id),]
p300.TypE <-  p300.gr.anno.enh[as.vector(idd),]




library(RColorBrewer)
library(dplyr)
pal <- brewer.pal(name = "YlGnBu", n = 9) %>% rev()
library(ggplot2)
library(ggthemes)

head(df.score)
###2d hist 
ggplot(df.score, aes(x = NB.score, y = CB.score)) +
    stat_binhex(bins = 100) +
    scale_fill_gradientn(colours = pal, trans = "log", breaks = c(1, 10, 100, 1000)) +
    theme_bw(base_size = 18) + facet_wrap(~chpSeq)
ggplot(df, aes(x = log2FC.exp, y = HiC.score)) +
    stat_binhex(bins = 100) + ylim(c(-2,4)) +
    scale_fill_gradientn(colours = pal, trans = "log", breaks = c(1, 10, 100, 1000)) +
    theme_bw(base_size = 18) + facet_wrap(~chpSeq)
###using tumor atlas


####BLACKLIST####

library(rtracklayer)
blacklist <- import.bed("encodeBlack.bed")

library(GenomicRanges)

black <- subsetByOverlaps(Atlas.idr, blacklist)

Atlas.idr.blackfilter <- import.bed("Atlas.idr.blackfilter.bed")

subsetByOverlaps(Atlas.idr.blackfilter, blacklist)




gel <- read.csv("../../bench/Ly7Ly1.x2.csv", header = T)


"tumorAtlas.bed"




###Tumors ####

atacT.load <- import.bed("tumorAtlas.filt.12.5kb.bed")
atacT.load <- summarizeOverlaps(atacT.load, samps$bamReads,singleEnd=FALSE, ignore.strand=TRUE, fragments=FALSE, mode="IntersectionNotEmpty")
save(atacT.load, file="atacT.load.Rdata")
