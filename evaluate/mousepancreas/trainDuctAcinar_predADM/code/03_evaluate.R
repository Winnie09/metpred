setwd('/home/whou10/data/whou/metpred/')
w <- readRDS('data/mousepancreas/wgbs/bs.rds')
source('/home/whou10/scratch16/whou10/resource/startup.R')
for (fn in list.files('evaluate/mousepancreas/trainDuctAcinar_predADM/res/')){
  print(fn)
  rdir1 = 'evaluate/mousepancreas/trainDuctAcinar_predADM/perf/pcc/'
  rdir2 = 'evaluate/mousepancreas/trainDuctAcinar_predADM/perf/mse/'
  dir.create(rdir1, recursive = T)
  dir.create(rdir2, recursive = T)
  
  pred <- readRDS(paste0('evaluate/mousepancreas/trainDuctAcinar_predADM/res/', fn))
  str(pred) 
  
  gs <- w[, colnames(pred)]  
  gs <- gs[rownames(pred), ] 
  
  str(gs)
  cc <- corfunc(pred, gs)  
  summary(cc)  
  saveRDS(cc, paste0(rdir1, '/', fn))
  
  mse <- sapply(1:nrow(pred), function(i){
    tmpv <- pred[i, ] - gs[i, ]
    mean(tmpv * tmpv)
  })
  names(mse) <- rownames(pred)
  saveRDS(mse, paste0(rdir2, '/', fn))
}

