library(matrixStats)
set.seed(12345)
m <- readRDS('/home-4/zji4@jhu.edu/scratch/encode_compiled/Nov19/wgbs/hg38mat/mat.rds')
m <- m/100
testid <- sample(colnames(m),16)
trainid <- setdiff(colnames(m),testid)
tm <- m[,testid]
m <- m[,trainid]
sdv <- rowSds(m)
id <- which(sdv > 0.3)
m <- m[id,]
tm <- tm[id,]
saveRDS(m,file='/home-4/zji4@jhu.edu/scratch/metpred/wgbs/data/wgbs/train.rds')
saveRDS(tm,file='/home-4/zji4@jhu.edu/scratch/metpred/wgbs/data/wgbs/test.rds')
