nd <- nrow(readRDS('/home-4/zji4@jhu.edu/scratch/metpred/data/data/mearray/beta.rds'))
m <- readRDS('/home-4/zji4@jhu.edu/scratch/metpred/data/data/expr/mat.rds')
m <- m[,c('A549','GM12878','IMR-90','K562')]
b <- readRDS('/home-4/zji4@jhu.edu/scratch/metpred/cellline/data/data/database/normexpr.rds')
intg <- intersect(row.names(m),row.names(b))
m <- m[intg,]
b <- b[intg,]
bm <- rowMeans(apply(b,2,sort))
rn <- row.names(m)
m <- apply(m,2,function(i) bm[rank(i)])
row.names(m) <- rn
clu <- readRDS('/home-4/zji4@jhu.edu/scratch/metpred/cellline/data/data/database/geneclu.rds')
ecl <- t(sapply(1:max(clu),function(i) {
  colMeans(m[names(clu)[clu==i],])
}))
set.seed(12345)
library(xgboost)
tecl <- t(ecl)

pred <- matrix(0,nrow=nd,ncol=ncol(m))
num <- ceiling(nd/10)
for (rid in 1:10) {
  print(rid)
  id <- intersect(1:nd,(num*(rid-1)+1):(num*rid))
  res <- readRDS(paste0('/home-4/zji4@jhu.edu/scratch/metpred/cellline/model/res/singlesite_logistic/',rid,'.rds'))
  pred[id,] <- t(sapply(res,predict,tecl))
}
pred[pred < 0] <- 0
pred[pred > 1] <- 1
saveRDS(pred,file='/home-4/zji4@jhu.edu/scratch/metpred/cellline/predict/singlesite/bulk/pred.rds')

