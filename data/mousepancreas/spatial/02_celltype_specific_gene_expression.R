library(Matrix)
library(nnls)
d <- as.matrix(readRDS('/home/whou10/data/whou/metpred/data/mousepancreas/spatial/normexpr.rds'))
w <- as.matrix(readRDS('/home/whou10/data/whou/metpred/data/mousepancreas/spatial/weight.rds'))

samp <- sub(':.*','',colnames(d))

## Within each sample,
## for each gene, find the non-negative linear least squares (NNLS)
## that solves min||celltype_weight * x - expr||_2. 
## Here, x is celltype-specific gene expression. 
m <- sapply(unique(samp),function(s) {
  id <- which(samp==s)
  sd <- d[,id]
  sw <- w[id,]
  cse <- t(apply(sd,1,function(sdr) nnls(sw,sdr)$x))
  colnames(cse) <- colnames(sw)
  cse
},simplify = F)

for (i in 1:length(m)) {
  colnames(m[[i]]) <- paste0(names(m)[i],':',colnames(m[[i]]))
}

m <- do.call(cbind,m)

saveRDS(m,file='/home/whou10/data/whou/metpred/data/mousepancreas/spatial/cse.rds')


