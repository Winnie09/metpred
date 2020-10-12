library(data.table)
m <- fread('/home-4/zji4@jhu.edu/scratch/scdata/GSE81861_singapore/geodownload/GSE81861_Cell_Line_COUNT.csv.gz',data.table = F)
row.names(m) <- sub('.*_','',m[,1])
m <- as.matrix(m[,-1])
load('/home-4/zji4@jhu.edu/scratch/resource/gn/res/hg19.rda')
mtg <- geneid[grepl('^MT-',genename)]
p <- colSums(m[intersect(row.names(m),mtg),])/colSums(m)
m <- m[,p < 0.2]
cell <- colnames(m)

m <- fread('/home-4/zji4@jhu.edu/scratch/scdata/GSE81861_singapore/geodownload/GSE81861_Cell_Line_FPKM.csv.gz',data.table = F)
row.names(m) <- sub('.*_','',m[,1])
m <- as.matrix(m[,-1])
m <- m[,cell]
b <- readRDS('/home-4/zji4@jhu.edu/scratch/metpred/cellline/data/data/database/geneclu.rds')
m <- m[intersect(row.names(m),names(b)),]
res <- log2(m + 1)
saveRDS(res,file='/home-4/zji4@jhu.edu/scratch/metpred/cellline/data/data/scexpr/log2FPKM.rds')
library(SAVER)
res <- saver(m,estimates.only=T,ncores=30)
res <- log2(res + 1)
saveRDS(res,file='/home-4/zji4@jhu.edu/scratch/metpred/cellline/data/data/scexpr/log2saver.rds')
#bm <- rowMeans(apply(b,2,sort))
#resr <- apply(res,2,function(i) bm[rank(i)])
#row.names(resr) <- row.names(res)

