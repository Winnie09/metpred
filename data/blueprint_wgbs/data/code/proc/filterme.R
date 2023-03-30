suppressMessages(library(GenomicRanges))
suppressMessages(library(BSgenome.Hsapiens.UCSC.hg38))
d <- readRDS('/home-4/zji4@jhu.edu/scratch/metpred/data/blueprint_wgbs/data/proc/me.rds')
r1 <- sub('_.*','',row.names(d))
r2 <- sub('.*_','',row.names(d))
id <- which(r1 %in% paste0('chr',c(1:22,'X')))
d <- d[id,]
r1 <- r1[id]
r2 <- r2[id]
r2 <- as.numeric(r2)
seq <- as.data.frame(getSeq(Hsapiens, r1, start = r2, end = r2 + 1))[,1]
id <- which(seq=='CG')
d <- d[id,]
saveRDS(d,'/home-4/zji4@jhu.edu/scratch/metpred/data/blueprint_wgbs/data/proc/filterme.rds')
