l <- readRDS('/home/whou10/data/whou/metpred/data/tcga/hg38/2024/array/mat/list.rds')
len <- sapply(l,length)
ulen <- unique(len)
for (ul in ulen) {
  print(ul)
  id <- which(len==ul)
  sl <- l[id]
  k <- sapply(2:length(sl),function(i) identical(names(sl[[1]]),names(sl[[i]])))
  nid <- which(!k)
  if (length(nid)==0) {
    print('pass')
  } else {
    for (i in nid+1) sl[[i]] <- sl[[i]][sort(names(sl[[i]]))]
    k <- sapply(2:length(sl),function(i) identical(names(sl[[1]]),names(sl[[i]])))
    if (sum(!k)==0) {print('pass')}
  }
  sl <- do.call(cbind,sl)
  k <- rowMeans(is.na(sl))
  sl <- sl[k < 1,]
  saveRDS(sl,file=paste0('/home/whou10/data/whou/metpred/data/tcga/hg38/2024/array/mat/sample_',nrow(sl),'.rds'))
}
