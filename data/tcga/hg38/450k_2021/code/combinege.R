ge <- readRDS('/home/whou10/data/whou/metpred/data/tcga/hg38/450k_2021/filelist/ge.rds')
me <- readRDS('/home/whou10/data/whou/metpred/data/tcga/hg38/450k_2021/filelist/me.rds')
int <- intersect(ge$patient,me$patient)
af <- list.files('/home/whou10/data/whou/metpred/data/tcga/hg38/450k_2021/ge',full.names = T)
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

saveRDS(fd,file='/home/whou10/data/whou/metpred/data/tcga/hg38/450k_2021/combine/ge.rds')

