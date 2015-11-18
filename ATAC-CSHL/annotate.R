library(ChIPpeakAnno)

#biocLite("org.Hs.eg.db")
library(org.Hs.eg.db)

setwd("TEACHING/")

atacpeaks <- toGRanges("ATAC_LY1_peaks.bed")

data(TSS.human.GRCh37)

anno <- annotatePeakInBatch(atacpeaks, AnnotationData=TSS.human.GRCh37, PeakLocForDistance="middle"
                            output="nearestLocation")

# add gene names
anno <- addGeneIDs(annotatedPeak=anno,
                   orgAnn="org.Hs.eg.db",
                   IDs2Add="symbol")

annoclose <- anno[abs(anno$distancetoFeature) < 2000]

unique(na.omit(as.vector(annoclose$symbol)))


#anno$peak[table(anno$peak)>1]
pie(table(anno$insideFeature))
sock <- IGVsocket()
#IGVgoto(sock,"chr1:247553693-247553876")
IGVgoto(sock,"chr1:117117370-117117553")
#> IGVload(sock, "ATAC-LY1.hg19.sorted.nodup.nonM.bam")


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