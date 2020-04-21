library(preprocessCore)
m <- readRDS('/home-4/zji4@jhu.edu/scratch/encode_compiled/Nov19/wgbs/hg38mat/mat.rds')
dm <- dimnames(m)
m <- normalize.quantiles(m)
dimnames(m) <- dm
saveRDS(m,file='/home-4/zji4@jhu.edu/scratch/metpred/data/wgbs/qn.rds')
