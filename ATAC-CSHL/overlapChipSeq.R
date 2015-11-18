library(ChIPpeakAnno)

#biocLite("org.Hs.eg.db")
library(org.Hs.eg.db)

##importing called peaks from macs2
extraCols_narrowPeak <- c(peakName = "character", DisScore = "numeric", strand="character", foldChange = "numeric", neg.logPv="numeric", neglogQval= "numeric", distPeakStart= "numeric")

LY1 <- import.bed("~/Projects/ATAC/sample_LY1/sample_LY1.narrow_peaks.narrowPeak", extraCol= extraCols_narrowPeak)

#Overlaps with chipSeq peaks
library(rtracklayer)

tracks <- read.csv("LY1.chipseq.csv", stringsAsFactors = F, header=TRUE)
setwd("chipSeq/")

chipSeqDB <- GRangesList()
for (i in seq(1:length(tracks$track))) {
    bedfile <- tracks$filename[i]
    bedname <- tracks$track[i]
    target <- tracks$target[i]
    gr <- import.bed(bedfile, trackLine = FALSE)
    names(gr) <- rep(target, length(gr))
    chipSeqDB[[target]] <- import.bed(bedfile, trackLine = FALSE)
    
}


#make a table of overlaps
olTable <- list()
for (i in seq(1:length((chipSeqDB)))) {
    xt <- ftable(!is.na(findOverlaps(unlist(chipSeqDB[i]), LY1, select="arbitrary")))
    #xtt <- addmargins(xt, 2, FUN = list(list(Total = sum)))
    tfname <- names(chipSeqDB)[i]
    xtt <- addmargins(xt, 2, FUN = list(list(Total = sum)))
    olTable[[tfname]] <- as.vector(xtt)
}


chip.df <- as.data.frame.list(olTable)
colnames(chip.df) <- c("Non-overlapping", "Overlapping", "Total")
chip.df$per.ol <- (chip.df$Overlapping / chip.df$Total) * 100
chip.df$chipExp <- row.names(chip.df)

chip.df.LY1 <- chip.df[order(-chip.df$per.ol),]
chip.df.LY1$chipExp <- reorder(chip.df.LY1$chipExp, chip.df.LY1$per.ol)
