cn <- commandArgs(trailingOnly = T)
d <- readRDS('/home-4/zji4@jhu.edu/scratch/metpred/data/data/mearray/beta.rds')
d <- d[,c('A549','GM12878','IMR-90','K562')]
pred <- readRDS(paste0('/home-4/zji4@jhu.edu/scratch/metpred/cellline/predict/singlesite/sc/pred_',cn,'.rds'))

sampcor <- sapply(1:ncol(d),function(i) cor(d[,i],pred[,i]))
sitecor <- sapply(1:nrow(d),function(i) cor(d[i,],pred[i,]))
sd <- apply(d,1,sd)
saveRDS(list(samplecor=sampcor,sitecor=sitecor,sitesd=sd),file=paste0('/home-4/zji4@jhu.edu/scratch/metpred/cellline/predict/singlesite/sc/cor_',cn,'.rds'))

