d <- readRDS('/home-4/zji4@jhu.edu/scratch/metpred/data/data/mearray/beta.rds')
ecl <- readRDS('/home-4/zji4@jhu.edu/scratch/metpred/cellline/data/data/database/clumean.rds')
d <- d[,colnames(ecl)]
sdv <- apply(d,1,sd)
id <- which(sdv > 0.2)
d <- d[id,]
saveRDS(d,file='/home-4/zji4@jhu.edu/scratch/metpred/cellline/data/data/metclu/submet.rds')
saveRDS(id,file='/home-4/zji4@jhu.edu/scratch/metpred/cellline/data/data/metclu/subid.rds')

d <- t(apply(d,1,scale))
clu <- kmeans(d,1000,iter.max=10000)$cluster
saveRDS(clu,file='/home-4/zji4@jhu.edu/scratch/metpred/cellline/data/data/metclu/clu1000.rds')
