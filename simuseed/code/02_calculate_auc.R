## for each seed, each setting, calculate all methods' perf (fdr,diff, auc)
rm(list=ls())
for (seed in 51:100){
  print(seed)
  source('/home-4/whou10@jhu.edu/work-zfs/whou10/trajectory_variability/function/01_function.R')
  ddir1 <- paste0('/home-4/whou10@jhu.edu/work-zfs/whou10/metpred/simuseed/res/', seed,'/')
  rdir1 <- paste0('/home-4/whou10@jhu.edu/work-zfs/whou10/metpred/simuseed/perf/', seed, '/sensfdr/')
  rdir2 <- paste0('/home-4/whou10@jhu.edu/work-zfs/whou10/metpred/simuseed/perf/', seed, '/perf/')
  af <- list.files(ddir1)
  #af <- af[!af %in% c('4_200_0.3_0.5')]
  dir.create(rdir2, recursive = TRUE)
  
  library(parallel)
  tmp <- lapply(af, function(f){
    if (!file.exists(paste0(rdir2, f, '.rds'))){
      print(f)
      dir.create(paste0(rdir1, f), recursive = TRUE)
      am <- list.files(paste0(ddir1, f))
      
      perf <- t(sapply(am, function(m){
        res <- readRDS(paste0(ddir1, f, '/', m))
        sensfdr <- SensFdr(grep('add:',rownames(res),value=T), res)
        saveRDS(sensfdr, paste0(rdir1, f, '/', m))
        AreaUnderSensFdr(sensfdr)
      }))
      rownames(perf) <- sub('_saver.*', '', rownames(perf))
      rownames(perf) <- sub('.rds.*', '', rownames(perf))
      saveRDS(perf, paste0(rdir2, f, '.rds'))
    }  
    return(0)
  })
}


