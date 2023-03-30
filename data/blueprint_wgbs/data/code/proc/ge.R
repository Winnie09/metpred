library(data.table)
af <- list.files('/home-4/zji4@jhu.edu/scratch/metpred/data/blueprint_wgbs/data/down/r')
d <- sapply(af,function(f) {
  tmp <- fread(paste0('/home-4/zji4@jhu.edu/scratch/metpred/data/blueprint_wgbs/data/down/r/',f),data.table = F)
  v <- tmp[,'TPM']
  names(v) <- sub('\\..*','',tmp[,1])
  v
})

saveRDS(d,file='/home-4/zji4@jhu.edu/scratch/metpred/data/blueprint_wgbs/data/proc/ge.rds')
