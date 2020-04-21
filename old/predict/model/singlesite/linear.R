library(parallel)
rid <- as.numeric(commandArgs(trailingOnly = T))
d <- readRDS('/home-4/zji4@jhu.edu/scratch/metpred/data/mearray/beta.rds')
d[d==0] <- min(d[d>0])
d <- log(d/(1-d))
d <- d[,readRDS('/home-4/zji4@jhu.edu/scratch/metpred/res/predict/label/trainid.rds')]

ecl <- readRDS('/home-4/zji4@jhu.edu/scratch/metpred/data/expr/clusterexpr.rds')
ecl <- ecl[,colnames(d)]

set.seed(12345)
library(xgboost)
tecl <- t(ecl)
num <- ceiling(nrow(d)/10)
id <- intersect(1:nrow(d),(num*(rid-1)+1):(num*rid))
mod <- lapply(id,function(i) {
  print(i)
  xgboost(tecl,d[i,],object='reg:linear',nrounds=10,verbose=F)
})
saveRDS(list(id,mod),file=paste0('/home-4/zji4@jhu.edu/scratch/metpred/res/predict/model/singlesite/linear/',rid,'.rds'))

