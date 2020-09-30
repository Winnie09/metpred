library(data.table)
d <- fread('/home-4/zji4@jhu.edu/scratch/metpred/hsmm/data/raw/GSE52529_fpkm_matrix.txt.gz',data.table = F)
row.names(d) <- d[,1]
d <- as.matrix(d[,-1])
saveRDS(d,file='/home-4/zji4@jhu.edu/scratch/metpred/hsmm/data/proc/fpkm.rds')
library(SAVER)
d <- log2(saver(d,estimates.only=T,ncores=20) + 1)
saveRDS(d,file='/home-4/zji4@jhu.edu/scratch/metpred/hsmm/data/proc/log2saver.rds')
