# setwd('/home-4/zji4@jhu.edu/scratch')
setwd('/scratch/users/whou10@jhu.edu/Wenpin')
profunc <- function(sgr,up,down) {
  start <- ifelse(strand(sgr)=='+',start(sgr)-up,end(sgr) + down)
  end <- ifelse(strand(sgr)=='+',start(sgr)-down,end(sgr) + up)
  GRanges(seqnames=as.character(seqnames(sgr)),IRanges(start=start,end=end),strand=strand(sgr))
}
suppressMessages(library(GenomicRanges))
suppressMessages(library(data.table))
gtf <- fread('/home-4/zji4@jhu.edu/scratch/resource/gtf/grch38.gtf',data.table = F)
gtf <- gtf[gtf[,3]=='gene',]
gid <- gsub('\"','',sub('gene_id ','',sapply(gtf[,9],function(i) strsplit(i,'; ')[[1]][1],USE.NAMES = F)))
gene <- GRanges(seqnames=gtf[,1],IRanges(start=gtf[,4],end=gtf[,5]),strand=gtf[,7])
names(gene) <- gid

gr <- readRDS('./encode_compiled/Nov19/wgbs/hg38filtermat/gr.rds') ## cpg site
m <- readRDS('./encode_compiled/Nov19/wgbs/hg38filtermat/mat.rds') ## cpg by celltype matrix
e <- readRDS('./metpred/data/data/expr/mat.rds') ## gene by celltype
e <- e[rowSums(e > 1) > 1,]
## match gene id
id = intersect(names(gene), rownames(e))
gene = gene[id]
e = e[id,]

pro <- profunc(gene[row.names(e),],up=500,down=0)
o <- as.matrix(findOverlaps(gr,pro))

lcv <- sapply(1:nrow(o),function(i) {
  print(i)
  cor(m[o[i,1],],e[o[i,2],])
})
saveRDS(list(m=m, e=e,o=o, lcv=lcv),'./metpred/plot/wgbs/plot/correlation.pd')

