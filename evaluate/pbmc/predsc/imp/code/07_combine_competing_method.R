af <- list.files('/home/whou10/data/whou/metpred/evaluate/pbmc/predsc/imp/competing/res/',pattern = 'permu_')
d <- sapply(af,function(f) {
  d <- readRDS(paste0('/home/whou10/data/whou/metpred/evaluate/pbmc/predsc/imp/competing/res/',f))
  cell <- sub(':.*','',colnames(d))
  cell[grepl('_t', cell)] <- 't_cells'
  sapply(sort(unique(cell)),function(i) {
    rowMeans(d[,cell==i])
  })
},simplify = F)
mean(sapply(2:length(d),function(i) identical(colnames(d[[i]]),colnames(d[[1]]))))
d <- do.call(rbind,d)
saveRDS(d,file='/home/whou10/data/whou/metpred/evaluate/pbmc/predsc/imp/competing/combine/permu.rds')


af <- list.files('/home/whou10/data/whou/metpred/evaluate/pbmc/predsc/imp/competing/res/',pattern = 'neargene_')
d <- sapply(af,function(f) {
  d <- readRDS(paste0('/home/whou10/data/whou/metpred/evaluate/pbmc/predsc/imp/competing/res/',f))
  cell <- sub(':.*','',colnames(d))
  cell[grepl('_t', cell)] <- 't_cells'
  sapply(sort(unique(cell)),function(i) {
    rowMeans(d[,cell==i])
  })
},simplify = F)
mean(sapply(2:length(d),function(i) identical(colnames(d[[i]]),colnames(d[[1]]))))
d <- do.call(rbind,d)
rownames(d) = sapply(rownames(d), function(i) sub(':','_',i), USE.NAMES = F)
saveRDS(d,file='/home/whou10/data/whou/metpred/evaluate/pbmc/predsc/imp/competing/combine/neargene.rds')


