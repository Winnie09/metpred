args <- commandArgs(trailingOnly = TRUE)
id <- as.numeric(args[1])
print('id ')
print(id)

alls <- list.files('/home/whou10/data/whou/metpred/evaluate/mousepancreas/pred/ramp/')
alls <- sub('.rds', '', alls) 
s <- alls[id]
print(s)
pred <- readRDS(paste0('/home/whou10/data/whou/metpred/evaluate/mousepancreas/pred/ramp/', s, '.rds'))

for (cutoff in seq(0.8, 0.95, 0.05)){
  print('cutoff ')
  print(cutoff)
  
  rdir1 <- paste0('/home/whou10/data/whou/metpred/evaluate/mousepancreas/pred/subset_weightSum_cutoff_spotList_', cutoff, '/')
  ##rdir2 <- paste0('/home/whou10/data/whou/metpred/evaluate/mousepancreas/goldstandard/', cutoff, '/')
  dir.create(rdir1, recursive = TRUE, showWarnings = F)
  ##dir.create(rdir2, recursive = TRUE)
  
  spot.list <- readRDS(paste0('/home/whou10/data/whou/metpred/evaluate/mousepancreas/goldstandard/res/weightSum_cutoff_', cutoff, '_spotList.rds'))
  print(s)
  spot <- spot.list[[s]]
  pred <- pred[, spot]
  saveRDS(pred, paste0(rdir1, s, '.rds'))
}



