ge <- readRDS('/home/whou10/data/whou/metpred/data/tcga/hg38/450k_2021/filelist/ge.rds')
me <- readRDS('/home/whou10/data/whou/metpred/data/tcga/hg38/450k_2021/filelist/me.rds')
int <- intersect(ge$patient,me$patient)
af <- list.files('/home/whou10/data/whou/metpred/data/tcga/hg38/450k_2021/me',full.names = T)
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
saveRDS(fd,file='/home/whou10/data/whou/metpred/data/tcga/hg38/450k_2021/combine/me.rds')



