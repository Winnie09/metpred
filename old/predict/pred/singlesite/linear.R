ecl <- readRDS('/home-4/zji4@jhu.edu/scratch/metpred/data/expr/clusterexpr.rds')
ecl <- ecl[,readRDS('/home-4/zji4@jhu.edu/scratch/metpred/res/predict/label/testid.rds')]
library(xgboost)
tecl <- t(ecl)
i <- as.numeric(commandArgs(trailingOnly = T))
res <- readRDS(paste0('/home-4/zji4@jhu.edu/scratch/metpred/res/predict/model/singlesite/linear/',i,'.rds'))
tmp <- t(sapply(1:length(res[[2]]),function(i) {print(i/length(res[[2]]));predict(res[[2]][[i]],tecl)}))
saveRDS(list(res[[1]],tmp),file=paste0('/home-4/zji4@jhu.edu/scratch/metpred/res/predict/pred/singlesite/linear/',i,'.rds'))
