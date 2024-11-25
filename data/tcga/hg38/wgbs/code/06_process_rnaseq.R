af <- list.files('/home/whou10/data/whou/metpred/data/tcga/hg38/wgbs/raw/download/rna',pattern = 'TCGA')
d <- sapply(af,function(f) {
  print(f)
  d <- data.table::fread(paste0('/home/whou10/data/whou/metpred/data/tcga/hg38/wgbs/raw/download/rna/',f),data.table=F)
  d <- d[nchar(d[,2]) > 0,]
  v <- d$tpm_unstranded
  names(v) <- d$gene_name
  v <- log2(v+1)
  v
},simplify = F)
unique(sapply(1:length(d),function(i) identical(names(d[[1]]),names(d[[i]]))))
d <- do.call(cbind,d)
saveRDS(d,file='/home/whou10/data/whou/metpred/data/tcga/hg38/wgbs/combine/ge.rds')
