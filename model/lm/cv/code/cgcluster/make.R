i <- as.numeric(commandArgs(trailingOnly = T))
load(paste0('/home-4/zji4@jhu.edu/scratch/metpred/model/lm/cv/res/lm_genecluster1000_geonly/',i,'.rda'))
d <- am
m <- rowMeans(d)
v <- (rowMeans(d*d) - m^2) * ((ncol(d) + 1)/ncol(d))
s <- sqrt(v)
d <- (d-m)/s
set.seed(12345)
clu <- kmeans(d,5000)$cluster
saveRDS(clu,file=paste0('/home-4/zji4@jhu.edu/scratch/metpred/model/lm/cv/res/cgcluster/',i,'.rds'))
