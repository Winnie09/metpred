library(parallel)
d <- readRDS('/home-4/zji4@jhu.edu/scratch/metpred/data/mearray/beta.rds')
d <- d[,readRDS('/home-4/zji4@jhu.edu/scratch/metpred/res/predict/label/trainid.rds')]

ecl <- readRDS('/home-4/zji4@jhu.edu/scratch/metpred/data/expr/clusterexpr.rds')
ecl <- ecl[,colnames(d)]

set.seed(12345)
library(xgboost)
tecl <- t(ecl)
mod <- mclapply(1:nrow(d),function(i) {
print(i)/nrow(d)
  xgboost(tecl,d[i,],object='reg:logistic',nrounds=10,verbose=F)
},mc.cores=detectCores()-1)
saveRDS(mod,file='/home-4/zji4@jhu.edu/scratch/metpred/res/predict/model/singlesite/logistic.rds')

