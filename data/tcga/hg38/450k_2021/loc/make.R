library(data.table)
library(rtracklayer)
suppressMessages(library(GenomicRanges))
d <- fread('/home/whou10/data/whou/metpred/data/tcga/hg38/450k_2021/loc/HumanMethylation450_15017482_v1-2.csv',data.table=F,skip=7,fill=T) ## https://support.illumina.com/downloads/humanmethylation450_15017482_v1-2_product_files.html
d <- d[d[,1]==d[,2] & d$Genome_Build=='37',]
d <- d[!is.na(d[,1]),]
g <- GRanges(seqnames=paste0('chr',d[,'CHR']),IRanges(start=d[,'MAPINFO'],end=d[,'MAPINFO']))
names(g) <- d[,1]
ch = import.chain('/home/whou10/data/whou/software/UCSCtools/liftover/hg19ToHg38.over.chain') ## https://hgdownload.soe.ucsc.edu/goldenpath/hg19/liftOver/
lf = liftOver(g, ch)
lf <- unlist(lf)

suppressMessages(library(BSgenome.Hsapiens.UCSC.hg38))
seq <- as.data.frame(getSeq(Hsapiens, as.character(seqnames(lf)), start = start(lf), end = start(lf) + 1))[,1]
id <- which(seq=='CG')
lf <- lf[id,]
v <- paste0(as.character(seqnames(lf)),'_',start(lf))
names(v) <- names(lf)
saveRDS(v,file='/home/whou10/data/whou/metpred/data/tcga/hg38/450k_2021/loc/hg38.rds')

