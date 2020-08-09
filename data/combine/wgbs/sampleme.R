d <- readRDS('/home-4/zji4@jhu.edu/scratch/metpred/data/combine/wgbs/filterme.rds')
cm <- rowMeans(d)
csv <- sqrt((rowMeans(d*d) - cm^2) / (ncol(d) - 1) * ncol(d))
d <- d[csv >= 0.2,]
set.seed(12345)
d <- d[sample(1:nrow(d),10000),]
saveRDS(d,file='/home-4/zji4@jhu.edu/scratch/metpred/data/combine/wgbs/sampleme.rds')
