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

gr <- readRDS('/home-4/zji4@jhu.edu/scratch/encode_compiled/Nov19/mearray/data/filtermat/gr.rds')
m <- readRDS('/home-4/zji4@jhu.edu/scratch/metpred/data/data/mearray/beta.rds')
e <- readRDS('/home-4/zji4@jhu.edu/scratch/metpred/data/data/expr/mat.rds')
e <- e[rowSums(e > 1) > 1,]
pro <- profunc(gene[row.names(e),],up=500,down=0)
o <- as.matrix(findOverlaps(gr,pro))

lcv <- sapply(1:nrow(o),function(i) {
  cor(m[o[i,1],],e[o[i,2],])
})
pdf('/home-4/zji4@jhu.edu/scratch/metpred/plot/plot/hist_acrossample.pdf')
hist(lcv)
dev.off()
scv <- sapply(1:ncol(m),function(i) cor(m[o[,1],i],e[o[,2],i]))
pdf('/home-4/zji4@jhu.edu/scratch/metpred/plot/plot/hist_acrossite.pdf')
hist(scv)
dev.off()

pdf('/home-4/zji4@jhu.edu/scratch/metpred/plot/plot/GM12878_acrosssite.pdf')
smoothScatter(m[o[,1],'GM12878']~e[o[,2],'GM12878'],xlab='Expression',ylab='Methylation')
dev.off()
pdf('/home-4/zji4@jhu.edu/scratch/metpred/plot/plot/acrossample.pdf')
plot(m[o[,1],][100000,]~e[o[,2],][100000,],xlab='Expression',ylab='Methylation')
dev.off()

