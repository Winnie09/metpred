timefunc <- function(i) {
  if (!grepl('-',i)) {
    t = as.numeric(times(i))
  } else {
    t = as.numeric(sub('-.*','',i))+as.numeric(times(sub('.*-','',i)))
  }
  round(60 * 24 * t,3)
}

setwd('/home/whou10/data/whou/metpred/evaluate/eff/out')
af = list.files(getwd())
af = af[3:4]
library(chron)
eff <- lapply(af,function(f){
  print(f)
  # d=readLines(paste0(f,'/out'))
  d=readLines(f)
  a = sapply(d[grep('whou10',d)+1],function(i) {
    tmp <- strsplit(i,' ')[[1]]
    tmp <- tmp[nchar(tmp) > 0]
    if ('COMPLETED' %in% tmp) {
      tmp[(length(tmp)-1):length(tmp)]  
    } 
    # else if ('FAILED' %in% tmp | 'RUNNING' %in% tmp) {
    #   c(NA,NA)
    # }
  },USE.NAMES = F)
  if (is.list(a)){
    a = a[!is.na(a)]
    a = do.call(cbind,a)
  } 
  a = a[,order(as.numeric(sub('K','',a[2,]))),drop=F]
  a[1,] <- sapply(a[1,],timefunc)
  a = cbind(sub('.txt','',f),a)
  colnames(a) = c('method',seq(100, (ncol(a)-1)*100, 100))
  a
})
names(eff) = af
alln = unique(unlist(sapply(eff, colnames))) 


eff = lapply(eff, function(i){
  tmp = matrix(NA, nrow = nrow(i), ncol = length(alln))
  dimnames(tmp) = list(rownames(i), alln)
  tmp[, colnames(i)] <- i
  tmp
})
res = do.call(rbind,eff)
time = res[seq(1,nrow(res),2),]
mem = res[seq(2,nrow(res),2),]
saveRDS(time,'/home/whou10/data/whou/metpred/evaluate/eff/res/time_include_allcpg.rds')
saveRDS(mem,'/home/whou10/data/whou/metpred/evaluate/eff/res/memory_include_allcpg.rds')


