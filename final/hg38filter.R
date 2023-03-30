d <- readRDS('/hpc/group/jilab/zj/met/final/fullwgbs_hg38.rds')
d <- d[!grepl('chrY',rownames(d)) & rowMeans(is.na(d)) < 1,]
d <- d[,colMeans(!is.na(d)) > 0.8]
samp <- sub('_.*','',colnames(d))
tab <- table(samp)
d <- d[,samp %in% names(tab)[tab >= 5]]
vaper <- rowSums(!is.na(d))
d <- d[vaper >= 50,]
saveRDS(d,file='/hpc/group/jilab/zj/met/final/wgbs_hg38.rds')

samp <- sub('_.*','',colnames(d))
setwd('/hpc/group/jilab/zj/met/combine/rna')
af <- list.files()
af <- sapply(unique(samp),function(i) grep(i,af,value=T))
m <- sapply(af,readRDS)
gene <- table(unlist(sapply(m,rownames)))
gene <- names(gene)[gene==length(m)]
m <- sapply(m,function(i) i[gene,])
m <- do.call(cbind,m)
m <- m[,colnames(d)]
m <- m[rowSums(m) > 0,]
saveRDS(m,file='/hpc/group/jilab/zj/met/final/rna.rds')

