rm(list=ls())
source('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/function/01_function.R')
ddir1 <- '/home-4/whou10@jhu.edu/work-zfs/whou10/metpred/simu/res/'
ddir2 <- '/home-4/whou10@jhu.edu/work-zfs/whou10/metpred/simu/data/data/'
rdir1 <- '/home-4/whou10@jhu.edu/work-zfs/whou10/metpred/simu/perf/sensfdr/'
rdir2 <- '/home-4/whou10@jhu.edu/work-zfs/whou10/metpred/simu/perf/perf/'
pdir <- '/home-4/whou10@jhu.edu/work-zfs/whou10/metpred/simu/plot/plot/'
pddir <- '/home-4/whou10@jhu.edu/work-zfs/whou10/metpred/simu/plot/pd/'

af <- list.files(rdir2)
#af <- af[!grepl('0.2.rds', af)]
perf <- t(sapply(af, function(f){
  tmp <- readRDS(paste0(rdir2, f))
  if (nrow(tmp) == 7){
    rownames(tmp) <- sub('.rds','',rownames(tmp))
    return(tmp[,1])
  } else {
    return(NULL)
  }
}))
if (is.list(perf)){
  id <- which(sapply(perf, length) != 0)
  perf <- perf[id]
  perf <- do.call(rbind, perf)
  rownames(perf) <- af[id]
}
  
colnames(perf) <- c('GLMM', 'limma', 'MAST', 'ourmethod', 'scDD', 't-test','wilcox')

library(reshape2)
pd1 <- melt(perf)
colnames(pd1) <- c('type', 'method', 'performance')
pd1 <- cbind(pd1, probRead = gsub('.*_','',sub('.rds','',pd1$type)))

perf <- t(sapply(af, function(f){
  tmp <- readRDS(paste0(rdir2, f))
  if (nrow(tmp) == 7){
    rownames(tmp) <- sub('.rds','',rownames(tmp))
    return(tmp[,2])
  } else {
    return(NULL)
  }
}))
if (is.list(perf)){
  id <- which(sapply(perf, length) != 0)
  perf <- perf[id]
  perf <- do.call(rbind, perf)
  rownames(perf) <- af[id]
}
colnames(perf) <- c('GLMM', 'limma', 'MAST', 'ourmethod', 'scDD', 't-test','wilcox')
library(reshape2)
pd2 <- melt(perf)
colnames(pd2) <- c('type', 'method', 'performance')
pd2 <- cbind(pd2, probRead = gsub('.*_','',sub('.rds','',pd2$type)))

saveRDS(pd1, paste0(pddir, 'perf_fdrdiff.rds')) ######### caution
saveRDS(pd2, paste0(pddir, 'perf_auc.rds')) ######### caution
