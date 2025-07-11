af <- list.files('/home/whou10/data/whou/metpred/data/tcga/hg38/2024/array/mat/',pattern = 'sample_')
for (f in af) {
  l <- readRDS(paste0('/home/whou10/data/whou/metpred/data/tcga/hg38/2024/array/mat/',f))
  m <- readRDS('/home/whou10/data/whou/metpred/data/tcga/hg38/2024/meta/combine/filtercase.rds')
  m <- m[[1]]
  rownames(m) <- m[,1]
  m <- m[colnames(l),]
  caseid <- m[,'Case ID']
  tab <- table(caseid)
  singleid <- names(tab)[tab==1]
  singlel <- l[,rownames(m)[m[,'Case ID'] %in% singleid]]
  colnames(singlel) <- m[colnames(singlel),'Case ID']
  
  multiid <- names(tab)[tab>1]
  multil <- l[,rownames(m)[m[,'Case ID'] %in% multiid]]
  k <- rowsum(t(multil),m[colnames(multil),'Case ID'],na.rm=T)
  tab <- table(m[colnames(multil),'Case ID'])
  k <- k/as.vector(tab[rownames(k)])
  print(identical(rownames(singlel),colnames(k)))
  d <- cbind(singlel,t(k))
  saveRDS(d,file=paste0('/home/whou10/data/whou/metpred/data/tcga/hg38/2024/array/mat/',sub('sample','case',f)))
  saveRDS(colnames(d),file=paste0('/home/whou10/data/whou/metpred/data/tcga/hg38/2024/array/mat/',sub('sample','casename',f)))
}


