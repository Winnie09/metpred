allf <- sub('.rds','',list.files('/home/whou10/data/whou/metpred/data/mousepancreas/spatial/procimpute/res/individual/'))
res <- lapply(allf,function(f) {
  print(f)
  sexpr <- readRDS(paste0("/home/whou10/data/whou/metpred/data/mousepancreas/spatial/procimpute/res/individual/",f,'.rds'))
})
d <- do.call(cbind, res)
saveRDS(d,file="/home/whou10/data/whou/metpred/data/mousepancreas/spatial/procimpute/res/all/imputednormexpr.rds")


