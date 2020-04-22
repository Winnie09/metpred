set.seed(123)
library(matrixStats)
# read in dna methylation data
d <- readRDS('/home-4/zji4@jhu.edu/scratch/metpred/data/data/mearray/beta.rds')
testid <- sample(colnames(d),10)
trainid <- setdiff(colnames(d),testid)
true <- d[,testid]
d <- d[,trainid]
# only retain > 0.2 sd CpG sites
id <- which(rowSds(d) > 0.3)
#id <- 1:nrow(d)
# trainrowid <- sample(id,10000)
# testrowid <- sample(setdiff(id,trainrowid),10000)
# testrowid <- trainrowid
# true <- true[testrowid,]
# d <- d[trainrowid,]
true <- true[id,]
d <- d[id,]
saveRDS(d,file='/home-4/zji4@jhu.edu/scratch/metpred/bulkcv//data/trainy.rds')
saveRDS(true,file='/home-4/zji4@jhu.edu/scratch/metpred/bulkcv//data/testy.rds')

suppressMessages(library(GenomicRanges))
suppressMessages(library(BSgenome.Hsapiens.UCSC.hg19))
gr <- readRDS('/home-4/zji4@jhu.edu/scratch/encode_compiled/Nov19/mearray/data/filtermat/gr.rds')
seq <- toupper(as.data.frame(getSeq(Hsapiens, as.character(seqnames(gr[id,])), start = start(gr[id,])-250, end = end(gr[id,]) + 251))[,1])
nchr <- nchar(seq[1])
chrmat <- do.call(rbind,strsplit(seq,''))
id <- apply(chrmat,2,function(i) length(unique(i)))
chrmat <- chrmat[,id > 1]
chrmat <- matrix(as.numeric(cbind(chrmat=='A',chrmat=='T',chrmat=='G',chrmat=='C')),nrow=nrow(chrmat))
saveRDS(chrmat,file='/home-4/zji4@jhu.edu/scratch/metpred/bulkcv//data/seqmat.rds')

m <- readRDS('/home-4/zji4@jhu.edu/scratch/metpred/data/data/expr/array.rds')
m <- m[,trainid]
m <- m[rowSums(m >= 5) >= 1,]
library(preprocessCore)
dn <- dimnames(m)
m <- normalize.quantiles(m)
dimnames(m) <- dn
b <- m
clu <- kmeans(t(apply(m,1,scale)),nrow(m)/100,iter.max=10000)$cluster
ecl <- t(sapply(1:max(clu),function(i) {
  colMeans(m[names(clu)[clu==i],])
}))

saveRDS(ecl,file='/home-4/zji4@jhu.edu/scratch/metpred/bulkcv//data/trainx.rds')

# cluster and average CpG sites
# scaled <- t(apply(d,1,scale))
# dclu <- kmeans(scaled,1000,iter.max = 10000)$cluster
# dmean <- t(sapply(1:max(dclu),function(i) {
#   colMeans(d[dclu==i,])
# }))

m <- readRDS('/home-4/zji4@jhu.edu/scratch/metpred/data/data/expr/array.rds')
m <- m[,testid]
intg <- intersect(row.names(m),row.names(b))
m <- m[intg,]
b <- b[intg,]
bm <- rowMeans(apply(b,2,sort))
rn <- row.names(m)
m <- apply(m,2,function(i) bm[rank(i)])
row.names(m) <- rn
clu <- clu[names(clu) %in% row.names(m)]
ecl <- t(sapply(1:max(clu),function(i) {
  colMeans(m[names(clu)[clu==i],])
}))
saveRDS(ecl,file='/home-4/zji4@jhu.edu/scratch/metpred/bulkcv//data/testx.rds')


