method = commandArgs(trailingOnly = T)[1] ## permu or neargene
print(method)

setwd('/home/whou10/data/whou/metpred/')
for (cvid in 1:5){
  print(cvid)
  for (fn in list.files(paste0('evaluate/mousepancreas/celltypespecific_cv/competing/', method, '/cv', cvid))){
    print(fn)
    rdir1 = paste0('evaluate/mousepancreas/celltypespecific_cv/perf/', method, '/pcc/cv', cvid)
    rdir2 = paste0('evaluate/mousepancreas/celltypespecific_cv/perf/', method, '/mse/cv', cvid)
    dir.create(rdir1, recursive = T)
    dir.create(rdir2, recursive = T)
    
    pred <- readRDS(paste0('evaluate/mousepancreas/celltypespecific_cv/competing/', method, '/cv', cvid, '/',fn))
    str(pred) ## num [1:10159500, 1:8]
    w <- readRDS('data/mousepancreas/wgbs/bs.rds')
    gs <- w[, colnames(pred)] # 
    gs <- gs[rownames(pred), ] # [1:955541, 1:8]
    source('/home/whou10/scratch16/whou10/resource/startup.R')
    str(gs) ## num [1:10159500, 1:8]
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
}
