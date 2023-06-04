s <- as.character(commandArgs(trailingOnly = T)[[1]])
print(s)

rdir <- '/home/whou10/data/whou/metpred/evaluate/mousepancreas/goldstandard/res/'
pdir <- '/home/whou10/data/whou/metpred/evaluate/mousepancreas/goldstandard/plot/'
w <- readRDS('/home/whou10/data/whou/metpred/data/mousepancreas/spatial/weight.rds')
str(w)

## generate gold standard, 
## for each spot, use weighted average wgbs measurements 
## predicted: each spots' predicted cpg-site-specific beta values
allf <- list.files('/home/whou10/data/whou/metpred/evaluate/mousepancreas/pred/ramp/')
allf <- sub('.rds','',allf)
bs <- readRDS('/home/whou10/data/whou/metpred/data/mousepancreas/wgbs/bs.rds')
alls <- gsub(':.*', '', colnames(bs))
cutoff <- 0.8
dir.create(paste0(rdir, 'goldstandard_weightedSum_', cutoff), recursive = T, showWarnings = F)

spot.list <- readRDS(paste0(rdir, 'weightSum_cutoff_', cutoff, '_spotList.rds'))
spot <- spot.list[[s]]

if (length(spot) > 0){
  dnam.s <- bs[, alls == s, drop = F]
  ct.s <- sub('.*:','', colnames(dnam.s))
  gs <- sapply(spot, function(spott){
    tmp <-  dnam.s %*% t(w[spott, ct.s, drop = FALSE])
    as.vector(tmp)
  })
  str(gs)  
  rownames(gs) <- rownames(dnam.s)
  saveRDS(gs, paste0(rdir, 'goldstandard_weightedSum_', cutoff, '/',  s, '.rds'))
}



