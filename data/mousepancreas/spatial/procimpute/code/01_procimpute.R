allf <- sub('.rds','',list.files('/home/whou10/data/whou/metpred/data/mousepancreas/spatial/imp/res/'))
res <- lapply(allf,function(f) {
  print(f)
  sexpr <- readRDS(paste0('/home/whou10/data/whou/metpred/data/mousepancreas/spatial/imp/res/',f,'.rds'))$estimate
  sexpr <- log2(sexpr + 1)
  saveRDS(sexpr,file=paste0("/home/whou10/data/whou/metpred/data/mousepancreas/spatial/procimpute/res/",f,'.rds'))
})

