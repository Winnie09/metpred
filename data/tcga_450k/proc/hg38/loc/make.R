library(data.table)
library(rtracklayer)
suppressMessages(library(GenomicRanges))
d <- fread('/home-4/zji4@jhu.edu/scratch/metpred/data/tcga_450k/proc/hg38/loc/HumanMethylation450_15017482_v1-2.csv',data.table=F,skip=7,fill=T)
d <- d[d[,1]==d[,2] & d$Genome_Build=='37',]
d <- d[!is.na(d[,1]),]
g <- GRanges(seqnames=paste0('chr',d[,'CHR']),IRanges(start=d[,'MAPINFO'],end=d[,'MAPINFO']))
names(g) <- d[,1]
ch = import.chain('/home-4/zji4@jhu.edu/scratch/software/UCSCtools/liftover/hg19ToHg38.over.chain')
lf = liftOver(g, ch)
lf <- unlist(lf)

suppressMessages(library(BSgenome.Hsapiens.UCSC.hg38))
seq <- as.data.frame(getSeq(Hsapiens, as.character(seqnames(lf)), start = start(lf), end = start(lf) + 1))[,1]
id <- which(seq=='CG')
lf <- lf[id,]
v <- paste0(as.character(seqnames(lf)),'_',start(lf))
names(v) <- names(lf)
saveRDS(v,file='/home-4/zji4@jhu.edu/scratch/metpred/data/tcga_450k/proc/hg38/loc/hg38.rds')
