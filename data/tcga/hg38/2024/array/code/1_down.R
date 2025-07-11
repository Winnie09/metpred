setwd('/home/whou10/data/whou/metpred/data/tcga/hg38/2024/array/raw')
d <- readRDS('/home/whou10/data/whou/metpred/data/tcga/hg38/2024/meta/combine/filtercase.rds')
d <- d[[1]]
for (f in d[,1])
  system(paste0('wget https://api.gdc.cancer.gov/data/',f))



