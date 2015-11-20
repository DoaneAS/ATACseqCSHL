# if needed
source("http://bioconductor.org/biocLite.R")
#biocLite("Gviz")
#biocLite("BiocUpgrade")
library(Rsamtools)
library(Gviz)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)
#options(ucscChromosomeNames=FALSE)
 
bamfile <- 'LY1_ATAC_chr18.bam'

# these are the coordinates on chr18 we want to visualize
start = 60790578-100000
end = 60986613+100000

# create a track based on BAM file 
bam_track <- DataTrack(range      = bamfile,
                       genome     = "hg19",
                       name       = "Coverage",
                       window     = -1,
                       chromosome = "chr18")

# plot tracks
plotTracks(bam_track, from = start, to = end, type="hist")

# add coordinate track
gtrack <- GenomeAxisTrack()

# add chromosome location track
itrack <- IdeogramTrack(genome = 'hg19', chromosome = 'chr18')

# plot all tracks
plotTracks(list(itrack, gtrack, bam_track), from = start, to = end, type="hist")


# now create a gene track based on UCSC knownGene annotation
library(TxDb.Hsapiens.UCSC.hg19.knownGene)

txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
grtrack <- GeneRegionTrack(txdb, genome = 'hg19', from = start, to = end,
                           transcriptAnnotation = "symbol",
                           geneSymbols=T,
                           chromosome = 'chr18', name = "UCSC known genes")

plotTracks(list(itrack, gtrack, bam_track, grtrack), from = start, to = end, type="hist")


# create a CpG  island track with data from UCSC
cpgIslands <- UcscTrack(genome='hg19', chromosome='chr18',
                         track="cpgIslandExt", from=start, to=end,
                         trackType="AnnotationTrack", start="chromStart",
                          end="chromEnd", id="name", shape="box",
                        fill="#006400", name="CGI")


#inlcude gene symbol

symbols <- unlist(mapIds(org.Hs.eg.db, gene(grtrack), "SYMBOL", "ENTREZID", multiVals = "first"))
symbol(grtrack) <- symbols[gene(grtrack)]

# plot all tracks
plotTracks(list(itrack, gtrack, bam_track, grtrack, cpgIslands), from = start, to = end, type="hist")

