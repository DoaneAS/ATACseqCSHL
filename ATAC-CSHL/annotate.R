library(ChIPpeakAnno)
source("https://bioconductor.org/biocLite.R")
biocLite("SRAdb")
biocLite("org.Hs.eg.db")
library(org.Hs.eg.db)
library(rtracklayer)
setwd("ATAC-CSHL")

##set up to read Bed 6+4 narrowPeak file as output by MACS2.1

extraCols_narrowPeak <- c(peakName = "character", DisScore = "numeric", strand="character", foldChange = "numeric", neg.logPv="numeric", neglogQval= "numeric", distPeakStart= "numeric")

atacpeaks <- import.bed("LY1.chr18.narrowPeak", extraCol= extraCols_narrowPeak)

library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb = TxDb.Hsapiens.UCSC.hg19.knownGene
ge <- genes(txdb, columns=c("tx_name", "gene_id", "tx_type"))

anno = annotatePeakInBatch(atacpeaks, FeatureLocForDistance = "TSS",  
                                     PeakLocForDistance="middle", select = "first", 
                                     output="nearestLocation",  AnnotationData=genes(txdb))

#add gene symbols ##
anno <- addGeneIDs(anno, silence=TRUE, orgAnn="org.Hs.eg.db", feature_id_type="entrez_id", IDs2Add=c("symbol"))


annoclose <- anno[abs(anno$distancetoFeature) < 2000]

unique(na.omit(as.vector(annoclose$symbol)))

### chromosomal regions ###

aCR<-assignChromosomeRegion(anno, nucleotideLevel=FALSE, 
                            proximal.promoter.cutoff=2000L,
                            immediate.downstream.cutoff=2000L,
                            precedence=c("Promoters", "immediateDownstream", 
                                         "fiveUTRs", "threeUTRs", 
                                         "Exons", "Introns"), 
                            TxDb=txdb)

#plot the percentages

barplot(aCR$percentage)


#### Annotate with chip-Seq in LY1 #####
### 

library(rtracklayer)

tracks <- read.table("LY1.chipseqAnno.txt", stringsAsFactors = F, header=TRUE)

chipSeqDB <- GRangesList()
for (i in seq(1:length(tracks$tf))) {
    bedfile <- tracks$fileName[i]
    bedname <- tracks$sampleName[i]
    target <- tracks$tf[i]
    gr <- import.bed(bedfile, trackLine = FALSE)
    names(gr) <- rep(target, length(gr))
    chipSeqDB[[target]] <- gr
    
}


##chr18


chipSeqDB18 <- endoapply(chipSeqDB, function(x) {
    unlist(split(x, seqnames(x))["chr18"])
})


#make a table of overlaps
olTable <- list()
for (i in seq(1:length((chipSeqDB18)))) {
    xt <- ftable(!is.na(findOverlaps(unlist(chipSeqDB18[i]), LY1, select="arbitrary")))
    #xtt <- addmargins(xt, 2, FUN = list(list(Total = sum)))
    tfname <- names(chipSeqDB18)[i]
    xtt <- addmargins(xt, 2, FUN = list(list(Total = sum)))
    olTable[[tfname]] <- as.vector(xtt)
}


chip.df <- as.data.frame.list(olTable)
colnames(chip.df) <- c("Non-overlapping", "Overlapping", "Total")
chip.df$per.ol <- (chip.df$Overlapping / chip.df$Total) * 100
chip.df$chipExp <- row.names(chip.df)

chip.df.LY1 <- chip.df[order(-chip.df$per.ol),]
chip.df.LY1$chipExp <- reorder(chip.df.LY1$chipExp, chip.df.LY1$per.ol)


chip.df.LY1$type <- NA
chip.df.LY1$type[grepl("H3k*", chip.df.LY1$chipExp)] <- "Marks"

chip.df.LY1$type[!grepl("H3k*", chip.df.LY1$chipExp)] <- "Transcription_Factors"

### plot the percentage overlaps ##
ggplot(chip.df.LY1, aes(y=per.ol, chipExp, fill=type)) + geom_bar(stat = "identity") + theme_few() + scale_fill_tableau("tableau10medium") + 
    coord_cartesian(ylim=c(0,100)) + theme(text=element_text(colour="grey20",size=18), 
                                           axis.text.x = element_text(angle=0,hjust=.5,vjust=-1,face="plain"),
                                           axis.text.y = element_text(angle=0,hjust=1,vjust=0,face="plain"),  
                                           axis.title.x = element_text(angle=0,hjust=.5,vjust=0,face="plain"),
                                           axis.title.y = element_text(angle=90,hjust=.5,vjust=.5,face="plain"))


######


















### IGV ####
#anno$peak[table(anno$peak)>1]
pie(table(anno$insideFeature))
sock <- IGVsocket()
#IGVgoto(sock,"chr1:247553693-247553876")
IGVgoto(sock,"chr1:117117370-117117553")
#> IGVload(sock, "ATAC-LY1.hg19.sorted.nodup.nonM.bam")
##

####chipSeq overlaps






ovly7 <- findOverlaps(atacpeaks, atacpeaksLY7)
atacpeaks_noly7 <- atacpeaks[-queryHits(ovly7)]

ovly1 <- findOverlaps(atacpeaksLY7, atacpeaks)
atacpeaksLY7_noly1 <- atacpeaksLY7[-queryHits(ovly1)]

library(rtracklayer)
export(resize(atacpeaksLY7[annoLY7only_close$peak], width=200, fix="center"), con="LY7only.bed", format="bed")

# annotate LY7 specific peaks

annoLY7only <- annotatePeakInBatch(atacpeaksLY7_noly1, AnnotationData=TSS.human.GRCh37,
                                   output="nearestLocation", maxgap=1000L)
annoLY7only <- addGeneIDs(annotatedPeak=annoLY7only,
                          orgAnn="org.Hs.eg.db",
                          IDs2Add="symbol")

annoLY7only[abs(annoLY7only$distancetoFeature) < 5000]$symbol

annoLY7only_close <- annoLY7only[abs(annoLY7only$distancetoFeature) < 5000]
#IRF4, HOXA7, SOX2

unique(na.omit(as.vector(annoLY7only_close$symbol)))