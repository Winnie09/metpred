af <- list.files('/home-4/zji4@jhu.edu/scratch/metpred/data/CEEHRC_wgbs/data/down/rna')
samp <- sub('.*\\.','',sub('\\.mRNA-Seq.*','',af))
ef <- sapply(1:92,function(i) {
  stf <- unique(samp)[i]
  file.exists(paste0('/home-4/zji4@jhu.edu/scratch/metpred/data/CEEHRC_wgbs/data/proc/rnasep/',stf,'.rds'))
}) 
cat(paste0('qsub make.sh ',which(!ef)),sep='\n')

