ge <- readRDS('/home-4/zji4@jhu.edu/scratch/metpred/data/tcga_450k/filelist/hg38/ge.rds')
me <- readRDS('/home-4/zji4@jhu.edu/scratch/metpred/data/tcga_450k/filelist/hg38/me.rds')
int <- intersect(ge$patient,me$patient)
af <- list.files('/home-4/zji4@jhu.edu/scratch/metpred/data/tcga_450k/proc/hg38/me',full.names = T)
d <- lapply(af,readRDS)
d <- do.call(cbind,d)
fd <- sapply(int,function(s) {
  tar <- me$file_id[me$patient==s]
  if (length(tar)==1) {
    d[,tar]
  } else {
    rowMeans(d[,tar])
  }
})
rn <- substr(row.names(fd),1,2)
fd <- fd[rn=='cg',]
fd <- fd[rowMeans(is.na(fd)) < 1,]
saveRDS(fd,file='/home-4/zji4@jhu.edu/scratch/metpred/data/tcga_450k/proc/hg38/combine/me.rds')


