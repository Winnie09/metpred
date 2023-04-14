d <- readRDS('/home/whou10/data/whou/metpred/data/encode/wgbs/filter/GRCh38.rds')
g <- readRDS('/home/whou10/data/whou/metpred/data/encode/rna/proc/GRCh38.rds')
k <- colMeans(is.na(d))
d <- d[,k < 0.5]
g <- g[,colnames(d)]
k <- rowMeans(is.na(d))
d <- d[k < 1,]
saveRDS(d,file='/home/whou10/data/whou/metpred/data/encode/final/GRCh38_wgbs.rds')
saveRDS(g,file='/home/whou10/data/whou/metpred/data/encode/final/GRCh38_rna.rds')
