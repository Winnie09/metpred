d <- readRDS('/home/whou10/data/whou/metpred/data/encode_wgbs/wgbs/proc/GRCh38.rds')
m <- readRDS('/home/whou10/data/whou/metpred/data/encode_wgbs/meta/combine/wgbs.rds')
m <- m[match(colnames(d),m$file),]
k <- sapply(unique(m$Biosample.summary),function(i) {
  id <- which(m$Biosample.summary==i)
  if (length(id)==1) {
    d[,id]
  } else {
    rowMeans(d[,id],na.rm=T)
  }
})
saveRDS(k,file='/home/whou10/data/whou/metpred/data/encode_wgbs/processed/GRCh38/wgbs.rds')


d <- readRDS('/home/whou10/data/whou/metpred/data/encode_wgbs/rna/proc/GRCh38.rds')
m <- readRDS('/home/whou10/data/whou/metpred/data/encode_wgbs/meta/combine/rna.rds')
m <- m[match(colnames(d),m$file),]
r <- sapply(unique(m$Biosample.summary),function(i) {
  id <- which(m$Biosample.summary==i)
  if (length(id)==1) {
    d[,id]
  } else {
    rowMeans(d[,id],na.rm=T)
  }
})

r <- r[,colnames(k)]

saveRDS(k,file='/home/whou10/data/whou/metpred/data/encode_wgbs/processed/GRCh38/wgbs.rds')
saveRDS(r,file='/home/whou10/data/whou/metpred/data/encode_wgbs/processed/GRCh38/rna.rds')

