# QC suing Shortreads package

# install
source("http://bioconductor.org/biocLite.R")
biocLite("ShortRead")
#biocLite("GenomicAlignments")

library(Rsamtools)
library(ShortRead)
library(lattice)
setwd("TEACHING")

# load reads
fq <- readFastq("reads.fastq.gz")

# look at data
fq 

# how many reads
length(fq)

# read #50000 
fq[50000]

# get sequence
sread(fq[50000])


quality(fq[50000])

# get quality scores
quals<-as(head(quality(fq[50000])), "matrix")

# look
quals

# translate into proba
10 ^ (-quals / 10)

# make boxplot
boxplot(as.matrix(PhredQuality(quality(fq))), outline=F, ylim=c(0,41), main="Per Cycle Read Quality", xlab="Cycle", ylab="Phred Quality", col=(c("gold","darkgreen")))

# plot nucleotide content
abc<-alphabetByCycle(sread(fq), alphabet=c("A","C","T","G"))
abc.df <- data.frame(Nucleotide=rownames(abc), Cycle=as.vector(col(abc)), Count=as.vector(abc))
print(xyplot(Count ~ Cycle, abc.df, group=Nucleotide, type="l", auto.key=list(lines=TRUE, points=FALSE))) 

# super fast way to get full report
qas <- qa("reads.fastq.gz")

# visualize it
browseURL(report(qas))

# trim reads from left using a sliding window (2 out 5 reads below "4" score)
fqt <- trimTailw(fq, 2, "4", halfwidth=2)

# check
fqt

# check new distribution of read length
hist(width(fqt))

# write trimmed reads to disk
writeFastq(fqt, "treads.fastq.gz", "w")


