suppressMessages(library(GenomicRanges))
suppressMessages(library(data.table))
gtf <- fread('/home-4/zji4@jhu.edu/scratch/resource/gtf/hg19.gtf',data.table = F)
gtf <- gtf[gtf[,3]=='gene',]
gid <- gsub('\"','',sub('gene_id ','',sapply(gtf[,9],function(i) strsplit(i,'; ')[[1]][1],USE.NAMES = F)))
gene <- GRanges(seqnames=gtf[,1],IRanges(start=gtf[,4],end=gtf[,5]),strand=gtf[,7])
names(gene) <- gid

gr <- readRDS('/home-4/zji4@jhu.edu/scratch/encode_compiled/Nov19/mearray/data/filtermat/gr.rds')
pro <- promoters(gene,upstream = 0,downstream = 1)

d <- distanceToNearest(gr,pro)
d <- cbind(as.matrix(d),unlist(mcols(d)))
#identical(1:nrow(d),unname(d[,1]))
#e <- readRDS('/home-4/zji4@jhu.edu/scratch/metpred/data/expr/mat.rds')
saveRDS(d,file='/home-4/zji4@jhu.edu/scratch/metpred/data/nearesttss/nearesttss.rds')
