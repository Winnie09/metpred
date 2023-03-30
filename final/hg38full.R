setwd('/hpc/group/jilab/zj/met/combine/wgbs')
d <- sapply(list.files(),readRDS)
gr <- readRDS('/hpc/group/jilab/zj/met/gr/hg19tohg38.rds')
gr <- gr[!gr[,2] %in% gr[duplicated(gr[,2]),2],]

for (i in grep('hg19',names(d),value=T)) {
  print(i)
  tmp <- d[[i]]
  tmp <- tmp[rownames(tmp) %in% gr[,1],]
  rownames(tmp) <- gr[match(rownames(tmp),gr[,1]),2]
  k <- setdiff(gr[,2],rownames(tmp))
  if (length(k) > 0) {
    m <- matrix(NA,nrow=length(k),ncol=ncol(tmp),dimnames = list(k,colnames(tmp)))
    d[[i]] <- rbind(tmp,m)[gr[,2],]  
  } else {
    d[[i]] <- tmp[gr[,2],]  
  }
  rm('m')
  rm('tmp')
  rm('k')
}

for (i in grep('hg38',names(d),value=T)) {
  print(i)
  tmp <- d[[i]]
  tmp <- tmp[rownames(tmp) %in% gr[,2],]
  k <- setdiff(gr[,2],rownames(tmp))
  if (length(k) > 0) {
    m <- matrix(NA,nrow=length(k),ncol=ncol(tmp),dimnames = list(k,colnames(tmp)))
    d[[i]] <- rbind(tmp,m)[gr[,2],]  
  } else {
    d[[i]] <- tmp[gr[,2],]  
  }
  rm('k')
  rm('m')
  rm('tmp')
}
d <- do.call(cbind,d)
colnames(d) <- sub('.*rds.','',colnames(d))

saveRDS(d,file='/hpc/group/jilab/zj/met/fullwgbs_hg38.rds')

