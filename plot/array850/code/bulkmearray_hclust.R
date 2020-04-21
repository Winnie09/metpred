library(ggplot2)
e <- readRDS('/home-4/zji4@jhu.edu/scratch/metpred/data/data/mearray/beta.rds')
e <- e[apply(e,1,sd) > 0.2,]
colnames(e) <- sub('_originated_from_.*','',sub('_female_.*','',sub('_male_.*','',colnames(e))))
dn <- dimnames(e)
e <- t(apply(e,1,scale))
dimnames(e) <- dn
hclu <- hclust(dist(t(e)))
pdf('/home-4/zji4@jhu.edu/scratch/metpred/plot/plot/bulkmearray_hclust.pdf',width=20,height=6)
plot(hclu)
dev.off()
