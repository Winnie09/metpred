library(rhdf5)
library(Matrix)
x = read_h5ad('GTEx_8_tissues_snRNAseq_atlas_071421.public_obs.h5ad')
d <- x$layers[['counts']]
m <- x$obs
identical(rownames(m),rownames(d))
tissue <- rowsum(as.matrix(d),m$tissue)
tissue <- log2(tissue/rowSums(tissue)*1e7+1)
tissue <- t(tissue)
saveRDS(tissue,file='tissue_pb.rds')

ct <- rowsum(as.matrix(d),paste0(m$tissue,':',m[,'Cell types level 2']))
ct <- log2(ct/rowSums(ct)*1e7+1)
ct <- t(ct)
saveRDS(ct,file='tissuect_pb.rds')

