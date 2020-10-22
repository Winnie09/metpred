rm(list=ls())
source('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/function/01_function.R')
ddir1 <- '/home-4/whou10@jhu.edu/work-zfs/whou10/metpred/simu/res/'
ddir2 <- '/home-4/whou10@jhu.edu/work-zfs/whou10/metpred/simu/data/data/'
rdir1 <- '/home-4/whou10@jhu.edu/work-zfs/whou10/metpred/simu/perf/sensfdr/'
rdir2 <- '/home-4/whou10@jhu.edu/work-zfs/whou10/metpred/simu/perf/perf/'
pdir <- '/home-4/whou10@jhu.edu/work-zfs/whou10/metpred/simu/plot/plot/'
pddir <- '/home-4/whou10@jhu.edu/work-zfs/whou10/metpred/simu/plot/pd/'
af <- list.files(ddir1)
#af <- af[!af %in% c('4_200_0.3_0.5')]
library(parallel)
tmp <- lapply(af, function(f){
  if (!file.exists(paste0(rdir2, f, '.rds'))){
    print(f)
    dir.create(paste0(rdir1, f), recursive = TRUE)
    am <- list.files(paste0(ddir1, f))
    perf <- t(sapply(am, function(m){
      res <- readRDS(paste0(ddir1, f, '/', m))
      sensfdr <- SensFdr(grep('add:',rownames(res),value=T), res)
      saveRDS(sensfdr, paste0(rdir1, f, '/', m, '.rds'))
      AreaUnderSensFdr(sensfdr)
    }))
    rownames(perf) <- sub('_saver.*', '', rownames(perf))
    saveRDS(perf, paste0(rdir2, f, '.rds'))
  }  
  return(0)
})

