profunc <- function(sgr,up,down) {
  start <- ifelse(strand(sgr)=='+',start(sgr)-up,end(sgr) + down)
  end <- ifelse(strand(sgr)=='+',start(sgr)-down,end(sgr) + up)
  GRanges(seqnames=as.character(seqnames(sgr)),IRanges(start=start,end=end),strand=strand(sgr))
}
suppressMessages(library(GenomicRanges))
suppressMessages(library(data.table))
gtf <- fread('/home-4/zji4@jhu.edu/scratch/resource/gtf/hg19.gtf',data.table = F)
gtf <- gtf[gtf[,3]=='gene',]
gid <- gsub('\"','',sub('gene_id ','',sapply(gtf[,9],function(i) strsplit(i,'; ')[[1]][1],USE.NAMES = F)))
gene <- GRanges(seqnames=gtf[,1],IRanges(start=gtf[,4],end=gtf[,5]),strand=gtf[,7])
names(gene) <- gid

gr <- readRDS('/home-4/zji4@jhu.edu/scratch/encode_compiled/Nov19/mearray/data/mat/gr.rds')
m <- readRDS('/home-4/zji4@jhu.edu/scratch/metpred/data/mearray/beta.rds')
e <- readRDS('/home-4/zji4@jhu.edu/scratch/metpred/data/expr/mat.rds')
e <- e[rowSums(e > 1) > 1,]
pro <- profunc(gene[row.names(e),],up=500,down=0)
o <- as.matrix(findOverlaps(gr,pro))

lcv <- sapply(1:nrow(o),function(i) {
  cor(m[o[i,1],],e[o[i,2],])
})

scv <- sapply(1:ncol(m),function(i) cor(m[o[,1],i],e[o[,2],i]))

