library(IlluminaHumanMethylationEPICanno.ilm10b3.hg19)
library(data.table)
library(rtracklayer)
suppressMessages(library(GenomicRanges))
data(Locations)
g <- as.data.frame(Locations)
g <- GRanges(seqnames=g[,1],IRanges(start=g[,2],end=g[,2]),strand=g[,3],names=rownames(g))
ch = import.chain('/home/whou10/data/whou/software/UCSCtools/liftover/hg19ToHg38.over.chain') ## https://hgdownload.soe.ucsc.edu/goldenpath/hg19/liftOver/
lf = liftOver(g, ch)
lf <- unlist(lf)

suppressMessages(library(BSgenome.Hsapiens.UCSC.hg38))
seq <- as.data.frame(getSeq(Hsapiens, as.character(seqnames(lf)), start = start(lf), end = start(lf) + 1))[,1]
id <- which(seq=='CG')
lf <- lf[id,]
v <- paste0(as.character(seqnames(lf)),'_',start(lf))
names(v) <- unname(unlist(mcols(lf)))
saveRDS(v,file='/home/whou10/data/whou/metpred/data/tcga/hg38/2024/array/pos/EPIC.rds')
