library(Matrix)
library(nnls)
alls <- sub('.rds','',list.files('/home/whou10/data/whou/metpred/data/mousepancreas/spatial/procimpute/res/individual/'))
for (s in alls){
  print(s)
  d <- readRDS(paste0('/home/whou10/data/whou/metpred/data/mousepancreas/spatial/procimpute/res/individual/', s, '.rds'))
  w <- as.matrix(readRDS('/home/whou10/data/whou/metpred/data/mousepancreas/spatial/proc/all/weight.rds'))
  w <- w[colnames(d), ]
  ## Within each sample,
  ## for each gene, find the non-negative linear least squares (NNLS)
  ## that solves min||celltype_weight * x - expr||_2. 
  ## Here, x is celltype-specific gene expression. 
  cse <- t(apply(d,1,function(sdr) nnls(w,sdr)$x))
  colnames(cse) <- colnames(w)
  colnames(cse) <- paste0(s,':',colnames(cse))
  saveRDS(cse,file=paste0('/home/whou10/data/whou/metpred/data/mousepancreas/spatial/imputed_nnls/res/individual/', s, '.rds'))
}

