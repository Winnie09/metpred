alls <- sub('.rds','',list.files('/home/whou10/data/whou/metpred/data/mousepancreas/spatial/imputed_nnls/res/individual/'))
d <- lapply(alls, function(s){
  tmp <- readRDS(paste0('/home/whou10/data/whou/metpred/data/mousepancreas/spatial/imputed_nnls/res/individual/', s, '.rds'))
})
d <- do.call(cbind, d)
saveRDS(d, '/home/whou10/data/whou/metpred/data/mousepancreas/spatial/imputed_nnls/res/all/cse.rds')

