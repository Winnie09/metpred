library(data.table)
af <- list.files('/home/whou10/data/whou/metpred/data/tcga/hg38/2024/ge/raw')
l <- sapply(af,function(f) {
  d <- fread(paste0('/home/whou10/data/whou/metpred/data/tcga/hg38/2024/ge/raw/',f),data.table=F)
  d <- d[d$gene_name!='',]
  v <- d$tpm_unstranded
  names(v) <- paste0(d$gene_name,':',d$gene_id)
  v
},simplify = F)

table(sapply(2:length(l),function(i) identical(names(l[[1]]),names(l[[i]]))))
l <- do.call(cbind,l)
saveRDS(l,file='/home/whou10/data/whou/metpred/data/tcga/hg38/2024/ge/mat/tpm_sample.rds')
