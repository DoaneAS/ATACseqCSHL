library(Rsamtools)

# set wd 
setwd("TEACHING")

bamfile <- 'LY1_ATAC_chr18.bam'

# optional
# sortBam(bamfile, "bam.sorted") # sort entries in bam file
# indexBam("bam.sorted.bam") # index reads in bam file

bam <- scanBam(bamfile)

# what did we get ?
class(bam)

# fields available
names(bam[[1]])

# let's look at 1 read 
bam[[1]]$seq[1]

# let's count nt frequency
alphabetFrequency(bam[[1]]$seq[1])

# lets's plot GC content
gcFunction = function(x) { 
  alf <- alphabetFrequency(x, as.prob=TRUE)
  rowSums(alf[,c("G", "C")])
}
readGC <- gcFunction(bam[[1]][["seq"]])
hist(readGC)

# MAPQ
bam[[1]]$mapq

# histogram of mapq
hist(bam[[1]]$mapq)

# EXERCISE: calculate % mapped reads with qual below 20?

# let's identify covered regions (peaks)

# read in 
library("GenomicAlignments")
aln <- readGAlignments('LY1_ATAC_chr18.bam')

# what is this ?
class(aln)

# get coverage
cover <- coverage(aln)

# look at chr18 coverage
cover$chr18

# identify islands of coverage >= 1 (view in IRange)
islands <- slice(cover, lower=1)

# identify max height in each island
islandPeakHeight <- viewMaxs(islands)

# identify width of each island
islandWidth <- width(islands)

# find peaks 
peaks <- islands[islandPeakHeight >= 2 & islandWidth >=100]

# look at them
peaks$chr18

# how many we have
length(peaks$chr18)
