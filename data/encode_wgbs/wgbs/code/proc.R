library(data.table)
meta <- readRDS('/home/whou10/data/whou/metpred/data/encode_wgbs/meta/combine/wgbs.rds')
conv <- c('GRCh38'='Homo sapiens','mm10'='Mus musculus')
aaf <- sub('.bed.gz','',list.files('/home/whou10/data/whou/metpred/data/encode_wgbs/wgbs/bed'))
mean(aaf %in% meta$file)
for (org in c('GRCh38','mm10')) {
  af <- intersect(aaf,meta[which(meta$Organism==conv[org]),'file'])
  cc = rep("NULL", 14)
  cc[c(1,6)] <- 'character'
  cc[c(2,11)] <- 'numeric'
  d <- sapply(af,function(f) {
    print(f)
    d <- fread(paste0('/home/whou10/data/whou/metpred/data/encode_wgbs/wgbs/bed/',f,'.bed.gz'),data.table=F,colClasses=cc)
    v <- d[,4]/100
    names(v) <- paste0(d[,1],':',ifelse(d[,3]=='+',d[,2]+1,d[,2]))
    v
  },simplify = F)
  u <- NULL
  for (id in 1:length(d)) u <- union(u,names(d[[id]]))
  print(length(u))
  m <- matrix(NA,nrow=length(u),ncol=length(d),dimnames=list(u,names(d)))
  for (n in names(d)) {
    m[names(d[[n]]),n] <- d[[n]]
  }
  saveRDS(m,file=paste0('/home/whou10/data/whou/metpred/data/encode_wgbs/wgbs/proc/',org,'.rds'))
}


