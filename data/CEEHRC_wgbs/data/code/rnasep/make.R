library(parallel)
suppressMessages(library(rtracklayer))
af <- list.files('/home-4/zji4@jhu.edu/scratch/metpred/data/CEEHRC_wgbs/data/down/rna')
samp <- sub('.*\\.','',sub('\\.mRNA-Seq.*','',af))
library(data.table)
gtf <- fread('/home-4/zji4@jhu.edu/scratch/resource/gtf/grch38.gtf',data.table = F)
gtf <- gtf[gtf[,3]=='gene',]
gn <- sub('\\..*','',sub('\".*','',sub('.*gene_id \"','',gtf[,9])))
gr <- GRanges(seqnames=gtf[,1],IRanges(start=gtf[,4],end=gtf[,5]))

stf <- unique(samp)[as.numeric(commandArgs(trailingOnly = T))]
tf <- af[samp==stf]
d1 <- import(paste0('/home-4/zji4@jhu.edu/scratch/metpred/data/CEEHRC_wgbs/data/down/rna/',grep('forward',tf,value=T)))
d2 <- import(paste0('/home-4/zji4@jhu.edu/scratch/metpred/data/CEEHRC_wgbs/data/down/rna/',grep('reverse',tf,value=T)))
v1 <- unlist(mcols(d1))
v2 <- unlist(mcols(d2))
print(c(mean(v1 > 0),mean(v2 < 0)))
d <- c(d1,d2)
v <- c(v1,-v2)
o <- as.matrix(findOverlaps(gr,d))

s <- tapply(v[o[,2]],list(o[,1]),mean)
names(s) <- gn[as.numeric(names(s))]

e <- log2(s + 1)
saveRDS(e,file=paste0('/home-4/zji4@jhu.edu/scratch/metpred/data/CEEHRC_wgbs/data/proc/rnasep/',stf,'.rds'))
