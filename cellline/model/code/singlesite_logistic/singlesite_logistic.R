library(parallel)
rid <- as.numeric(commandArgs(trailingOnly = T))
d <- readRDS('/home-4/zji4@jhu.edu/scratch/metpred/data/data/mearray/beta.rds')
ecl <- readRDS('/home-4/zji4@jhu.edu/scratch/metpred/cellline/data/data/database/clumean.rds')
d <- d[,colnames(ecl)]

set.seed(12345)
library(xgboost)
tecl <- t(ecl)
num <- ceiling(nrow(d)/10)
id <- intersect(1:nrow(d),(num*(rid-1)+1):(num*rid))
mod <- lapply(id,function(i) {
print(i)
  xgboost(tecl,d[i,],object='reg:logistic',nrounds=10,verbose=F)
})
saveRDS(mod,file=paste0('/home-4/zji4@jhu.edu/scratch/metpred/cellline/model/res/singlesite_logistic/',rid,'.rds'))

