library(GenomicRanges) 
library(Rsamtools) 
library(rtracklayer) 
library(GenomicAlignments)

setwd("TEACHING")

anno

# set up BAM scanning parameters
params <- ScanBamParam(flag=scanBamFlag(isUnmappedQuery=FALSE))
  
# scan BAM file
bam <- readGAlignments("LY1_ATAC_chr18.bam", param=params)

# calculate counts
counts <- summarizeOverlaps( features=anno,
                             reads=bam,
                             mode="Union",
                             ignore.strand=TRUE,
                             SingleEnd=FALSE,
                             param=params) 

# what did we get ?
counts

# make count table
count_table <- assays(counts, withDimnames=TRUE)$counts

# EXERCISE 1: how many counts for BCL2 ?
bcl2 <- counts[mcols(counts)$symbol %in% "BCL2",]
assay(bcl2)

# EXERCISE 2: convert to count per million reads (cpm)

# EXERCISE 3: calculate RPKM ?

# write count table to disk
write.table(count_table,
            file="counts_table.txt", sep = "\t", 
row.names=TRUE, col.names=FALSE, quote=F)

1000000*1000*count_table/sum(width(reduce(genes)))/sum(count_table)
