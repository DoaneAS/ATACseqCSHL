source("https://bioconductor.org/biocLite.R")
biocLite()
biocLite("ChIPpeakAnno")
library(ChIPpeakAnno)
biocLite("org.Hs.eg.db")
biocLite("TxDb.Hsapiens.UCSC.hg19.knownGene")
library(org.Hs.eg.db)
library(rtracklayer)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)

library(ggplot2)

source("atacseqFuncitons.R")

setwd("ATAC-CSHL")

##set up to read Bed 6+4 narrowPeak file as output by MACS2.1

extraCols_narrowPeak <- c(peakName = "character", DisScore = "numeric", strand="character", foldChange = "numeric", neg.logPv="numeric", neglogQval= "numeric", distPeakStart= "numeric")

#import the file
atacpeaks <- import.bed("LY1.chr18.narrowPeak", extraCol= extraCols_narrowPeak)

#load genome annotation
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb = TxDb.Hsapiens.UCSC.hg19.knownGene
ge <- genes(txdb, columns=c("tx_name", "gene_id", "tx_type"))


## assign each peak to a gene according to the nearest TSS
#what are implications of doing this?  when would our assignments tend to be more or less accurate?

anno = annotatePeakInBatch(atacpeaks, FeatureLocForDistance = "TSS",  
                                     PeakLocForDistance="middle", select = "first", 
                                     output="nearestLocation",  AnnotationData=genes(txdb))

#add gene symbols ##

anno <- addGeneIDs(anno, silence=TRUE, orgAnn="org.Hs.eg.db", feature_id_type="entrez_id", IDs2Add=c("symbol"))


#what are the genes that are assigned to peaks and located with <= 2kb
annoclose <- anno[abs(anno$distancetoFeature) < 2000]

unique(na.omit(as.vector(annoclose$symbol)))

#characterize the overlap with respect to genes
pie(table(anno$insideFeature))

### chromosomal regions ###

aCR<-assignChromosomeRegion(anno, nucleotideLevel=FALSE, 
                            proximal.promoter.cutoff=2500L,
                            immediate.downstream.cutoff=2500L,
                            precedence=c("Promoters", "immediateDownstream", 
                                         "fiveUTRs", "threeUTRs", 
                                         "Exons", "Introns"), 
                            TxDb=txdb)

#plot the percentages

barplot(aCR$percentage)



##################################################################################
####    HOW WELL DO ATAC_SEQ PEAKS CAPTURE REGIONS IDENTIFIED BY CHIP-SEQ    #####
#### Annotate with chip-Seq peaks in LY1 from GSE 29282                      #####
##################################################################################



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

## limit to peaks on chr18 only ##


chipSeqDB18 <- endoapply(chipSeqDB, function(x) {
    unlist(split(x, seqnames(x))["chr18"])
})


#make a table of overlaps
olTable <- list()
for (i in seq(1:length((chipSeqDB18)))) {
    a <- ftable(!is.na(findOverlaps(unlist(chipSeqDB18[i]), anno, select="arbitrary")))
    tfname <- names(chipSeqDB18)[i]
    b <- addmargins(a, 2, FUN = list(list(Total = sum)))
    olTable[[tfname]] <- as.vector(b)
}

#make a data frame of the overlaps 
chip.df <- as.data.frame.list(olTable)
#name the columns
colnames(chip.df) <- c("Non-overlapping", "Overlapping", "Total")

#compute percentage overlap
chip.df$per.ol <- (chip.df$Overlapping / chip.df$Total) * 100
chip.df$chipExp <- row.names(chip.df)

#order dataframe by precent overlap
chip.df.LY1 <- chip.df[order(-chip.df$per.ol),]

# order factor levels by percent overlap
chip.df.LY1$chipExp <- reorder(chip.df.LY1$chipExp, chip.df.LY1$per.ol)

# make a categorical factor to distinguish histone marjs and transcription factors ##
chip.df.LY1$type <- NA
chip.df.LY1$type[grepl("H3k*", chip.df.LY1$chipExp)] <- "Marks"

chip.df.LY1$type[!grepl("H3k*", chip.df.LY1$chipExp)] <- "Transcription_Factors"

### plot the percentage overlaps ##
library(ggplot2)

(ggplot(chip.df.LY1, aes(y=per.ol, chipExp, fill=type)) + geom_bar(stat = "identity") + theme_minimal() + scale_fill_brewer(palette ="Accent") + 
    coord_cartesian(ylim=c(0,100)) + theme(text=element_text(colour="grey20",size=18), 
                                           axis.text.x = element_text(angle=0,hjust=.5,vjust=-1,face="plain"),
                                           axis.text.y = element_text(angle=0,hjust=1,vjust=0,face="plain"),  
                                           axis.title.x = element_text(angle=0,hjust=.5,vjust=0,face="plain"),
                                           axis.title.y = element_text(angle=90,hjust=.5,vjust=.5,face="plain")))

###
####  To GVIZ.R #################














