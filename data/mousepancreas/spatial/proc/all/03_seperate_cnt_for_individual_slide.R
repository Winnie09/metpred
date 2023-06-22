cnt <- readRDS(file='/home/whou10/data/whou/metpred/data/mousepancreas/spatial/proc/all/count.rds')
samp <- sub(':.*', '', colnames(cnt))
for (s in unique(samp)){
  cs <- as.matrix(cnt[, samp == s])
  saveRDS(cs, paste0('/home/whou10/data/whou/metpred/data/mousepancreas/spatial/proc/individual/cnt/', s, '.rds'))
}


