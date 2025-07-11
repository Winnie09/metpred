library(data.table)
af <- list.files('/home/whou10/data/whou/metpred/data/tcga/hg38/2024/array/raw')
l <- sapply(af,function(f) {
  d <- fread(paste0('/home/whou10/data/whou/metpred/data/tcga/hg38/2024/array/raw/',f),data.table=F)
  v <- d[,2]
  names(v) <- d[,1]
  v[grepl('^cg',names(v))]
},simplify = F)
saveRDS(l,file='/home/whou10/data/whou/metpred/data/tcga/hg38/2024/array/mat/list.rds')

