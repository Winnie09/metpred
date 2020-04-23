m <- readRDS('/home-4/zji4@jhu.edu/scratch/metpred/data/data/expr/mat.rds')
m <- m[,!colnames(m) %in% c('A549','GM12878','IMR-90','K562')]
m <- m[rowSums(m >= 1) >= 1,]
library(preprocessCore)
dn <- dimnames(m)
m <- normalize.quantiles(m)
dimnames(m) <- dn
saveRDS(m,file='/home-4/zji4@jhu.edu/scratch/metpred/cellline/data/data/database/normexpr.rds')
clu <- kmeans(t(apply(m,1,scale)),nrow(m)/100,iter.max=10000)$cluster
saveRDS(clu,file='/home-4/zji4@jhu.edu/scratch/metpred/cellline/data/data/database/geneclu.rds')
ecl <- t(sapply(1:max(clu),function(i) {
  colMeans(m[names(clu)[clu==i],])
}))
saveRDS(ecl,file='/home-4/zji4@jhu.edu/scratch/metpred/cellline/data/data/database/clumean.rds')
