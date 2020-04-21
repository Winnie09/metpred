library(parallel)
cellnum <- as.numeric(commandArgs(trailingOnly = T))
nd <- nrow(readRDS('/home-4/zji4@jhu.edu/scratch/metpred/cellline/data/data/metclu/submet.rds'))
m <- readRDS('/home-4/zji4@jhu.edu/scratch/metpred/cellline/data/data/scexpr/log2saver.rds')
ct <- sapply(colnames(m),function(i) strsplit(i,'__')[[1]][2],USE.NAMES = F)
ct <- sub('_.*','',ct)
ct[ct=='IMR90'] <- 'IMR-90'
m <- sapply(c('A549','GM12878','IMR-90','K562'),function(sct) {
  rowMeans(m[,sample(which(ct==sct),cellnum),drop=F])
})
b <- readRDS('/home-4/zji4@jhu.edu/scratch/metpred/cellline/data/data/database/normexpr.rds')
intg <- intersect(row.names(m),row.names(b))
m <- m[intg,]
b <- b[intg,]
bm <- rowMeans(apply(b,2,sort))
rn <- row.names(m)
m <- apply(m,2,function(i) bm[rank(i)])
row.names(m) <- rn
clu <- readRDS('/home-4/zji4@jhu.edu/scratch/metpred/cellline/data/data/database/geneclu.rds')
clu <- clu[names(clu) %in% row.names(m)]
ecl <- t(sapply(1:max(clu),function(i) {
  colMeans(m[names(clu)[clu==i],])
}))

set.seed(12345)
library(xgboost)
tecl <- t(ecl)

pred <- matrix(0,nrow=nd,ncol=ncol(m))
res <- readRDS(paste0('/home-4/zji4@jhu.edu/scratch/metpred/cellline/model/res/combinesite_logistic_tune/comb.rds'))
tmp <- mclapply(1:length(res),function(i) predict(res[[i]],tecl),mc.cores=20)
dclu <- readRDS('/home-4/zji4@jhu.edu/scratch/metpred/cellline/data/data/metclu/clu1000.rds')
tmp <- do.call(rbind,tmp)
pred <- tmp[dclu,]

pred[pred < 0] <- 0
pred[pred > 1] <- 1
saveRDS(pred,file=paste0('/home-4/zji4@jhu.edu/scratch/metpred/cellline/predict/singlesite/sc/predclu_',cellnum,'.rds'))


