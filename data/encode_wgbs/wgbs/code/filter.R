library(BSgenome.Hsapiens.UCSC.hg38)
library(GenomicRanges)
d <- readRDS('/home/whou10/data/whou/metpred/data/encode/wgbs/proc/GRCh38.rds')
s <- sub(':.*','',rownames(d))
d <- d[nchar(s) <= 5 & !s %in% c('chrY','chrM'),]
s <- sub(':.*','',rownames(d))
l <- as.numeric(sub('.*:','',rownames(d)))
g <- GRanges(seqnames=s,IRanges(start=l,end=l+1))
s <- as.character(Views(Hsapiens, g))

meta <- readRDS('/home/whou10/data/whou/metpred/data/encode/meta/combine/wgbs.rds')

type <- meta[match(colnames(d),meta[,1]),3]
md <- sapply(unique(type),function(i) {
  id <- which(type==i)
  if (length(id)==1) {
    d[,id]    
  } else {
    rowMeans(d[,id],na.rm=T)  
  }
})
saveRDS(md,file='/home/whou10/data/whou/metpred/data/encode/wgbs/filter/GRCh38.rds')



library(BSgenome.Mmusculus.UCSC.mm10)
library(GenomicRanges)
d <- readRDS('/home/whou10/data/whou/metpred/data/encode/wgbs/proc/mm10.rds')
s <- sub(':.*','',rownames(d))
d <- d[nchar(s) <= 5 & !s %in% c('chrY','chrM'),]
s <- sub(':.*','',rownames(d))
l <- as.numeric(sub('.*:','',rownames(d)))
g <- GRanges(seqnames=s,IRanges(start=l,end=l+1))
s <- as.character(Views(Mmusculus, g))

meta <- readRDS('/home/whou10/data/whou/metpred/data/encode/meta/combine/wgbs.rds')

type <- meta[match(colnames(d),meta[,1]),3]
colnames(d) <- type
saveRDS(d,file='/home/whou10/data/whou/metpred/data/encode/wgbs/filter/mm10.rds')

