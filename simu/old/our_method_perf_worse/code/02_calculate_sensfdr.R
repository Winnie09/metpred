rm(list=ls())
source('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/function/01_function.R')
ddir1 <- '/home-4/whou10@jhu.edu/work-zfs/whou10/metpred/simu/res/'
ddir2 <- '/home-4/whou10@jhu.edu/work-zfs/whou10/metpred/simu/data/data/'
rdir1 <- '/home-4/whou10@jhu.edu/work-zfs/whou10/metpred/simu/perf/sensfdr/'
rdir2 <- '/home-4/whou10@jhu.edu/work-zfs/whou10/metpred/simu/perf/perf/'
pdir <- '/home-4/whou10@jhu.edu/work-zfs/whou10/metpred/simu/plot/plot/'
pddir <- '/home-4/whou10@jhu.edu/work-zfs/whou10/metpred/simu/plot/pd/'
af <- list.files(ddir1)

library(parallel)
tmp <- lapply(af, function(f){
  if (!file.exists(paste0(rdir2, f, '.rds'))){
    print(f)
    dir.create(paste0(rdir1, f))
    am <- list.files(paste0(ddir1, f))
    perf <- t(sapply(am, function(m){
      res <- readRDS(paste0(ddir1, f, '/', m))
      sensfdr <- SensFdr(grep('add:',rownames(res),value=T), res)
      saveRDS(sensfdr, paste0(rdir1, f, '/', m, '.rds'))
      AreaUnderSensFdr(sensfdr)
    }))
    saveRDS(perf, paste0(rdir2, f, '.rds'))
  }  
  return(0)
})

