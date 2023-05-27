library(data.table)
d <- fread('GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.gz',data.table=F)
gn <- paste0(sub('\\..*','',d[,1]),':',d[,2])
d <- as.matrix(d[,-c(1:2)])
rownames(d) <- gn
d <- log2(d+1)
saveRDS(d,file='log2tpm.rds')
