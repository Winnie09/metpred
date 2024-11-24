ge <- readRDS('/home-4/zji4@jhu.edu/scratch/metpred/data/tcga_450k/filelist/hg38/ge.rds')
me <- readRDS('/home-4/zji4@jhu.edu/scratch/metpred/data/tcga_450k/filelist/hg38/me.rds')
int <- intersect(ge$patient,me$patient)
af <- list.files('/home-4/zji4@jhu.edu/scratch/metpred/data/tcga_450k/proc/hg38/ge',full.names = T)
d <- lapply(af,readRDS)
d <- do.call(cbind,d)
d <-log2(d+1)
fd <- sapply(int,function(s) {
  tar <- ge$file_id[ge$patient==s]
  if (length(tar)==1) {
    d[,tar]
  } else {
    rowMeans(d[,tar])
  }
})
fd <- fd[rowMeans(fd) > 0,]
saveRDS(fd,file='/home-4/zji4@jhu.edu/scratch/metpred/data/tcga_450k/proc/hg38/combine/ge.rds')

