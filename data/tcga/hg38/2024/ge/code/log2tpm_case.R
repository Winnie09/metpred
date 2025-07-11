l <- readRDS('/home/whou10/data/whou/metpred/data/tcga/hg38/2024/ge/mat/tpm_sample.rds')
m <- readRDS('/home/whou10/data/whou/metpred/data/tcga/hg38/2024/meta/combine/filtercase.rds')
m <- m[[2]]
rownames(m) <- m[,1]
k <- rowsum(t(l),m[colnames(l),'Case ID'])
tab <- table(m[colnames(l),'Case ID'])
k <- k/as.vector(tab[rownames(k)])
k <- log2(t(k) + 1)
saveRDS(k,file='/home/whou10/data/whou/metpred/data/tcga/hg38/2024/ge/mat/log2tpm_case.rds')
