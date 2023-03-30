library(data.table)
af <- list.files('/home-4/zji4@jhu.edu/scratch/metpred/data/GSE108066/data/raw/wgbs')
names(af) <- paste0('GSE108066_',sub('GSM[0-9]*_','',sub('_CpG_WGBS.txt.gz','',af)))
cn <- colnames(readRDS('/home-4/zji4@jhu.edu/scratch/metpred/data/GSE108066/data/proc/rna.rds'))
af <- af[cn]
m <- sapply(af,function(f) {
  d <- fread(paste0('/home-4/zji4@jhu.edu/scratch/metpred/data/GSE108066/data/raw/wgbs/',f),data.table=F)
  n <- paste0(d[,1],'_',d[,2])
  v <- d[,3]/d[,4]
  names(v) <- n
  v
})
colnames(m) <- names(af)
saveRDS(m,file='/home-4/zji4@jhu.edu/scratch/metpred/data/GSE108066/data/proc/wgbs.rds')
