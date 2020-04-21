library(ggplot2)
e <- readRDS('/home-4/zji4@jhu.edu/scratch/metpred/data/data/expr/mat.rds')
e <- e[rowSums(e >= 1) >= 1,]
library(preprocessCore)
dn <- dimnames(e)
e <- normalize.quantiles(e)
dimnames(e) <- dn
e <- e[apply(e,1,sd)/rowMeans(e) > 0.5,]
colnames(e) <- sub('_originated_from_.*','',sub('_female_.*','',sub('_male_.*','',colnames(e))))
dn <- dimnames(e)
e <- t(apply(e,1,scale))
dimnames(e) <- dn
hclu <- hclust(dist(t(e)))
pdf('/home-4/zji4@jhu.edu/scratch/metpred/plot/plot/bulkrnaseq_hclust.pdf',width=20,height=6)
plot(hclu)
dev.off()
