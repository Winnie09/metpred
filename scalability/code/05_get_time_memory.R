timefunc <- function(i) {
  if (!grepl('-',i)) {
    t = as.numeric(times(i))
  } else {
    t = as.numeric(sub('-.*','',i))+as.numeric(times(sub('.*-','',i)))
  }
  round(60 * 24 * t,3)
}

setwd('/home-4/whou10@jhu.edu/work-zfs/whou10/metpred/scalability/code/comp/')
allmtd = list.files(getwd())
mtd = allmtd[1]

library(chron)
res <- lapply(allmtd, function(mtd){
  ddir <- paste0(getwd(), '/', mtd, '/')
  af <- list.files(ddir)
  
  eff <- lapply(af,function(f){
    print(f)
    d=readLines(paste0(ddir, f,'/out'))
    a = sapply(d[grep('whou10',d)+1],function(i) {
      tmp <- strsplit(i,' ')[[1]]
      tmp <- tmp[nchar(tmp) > 0]
      if ('COMPLETED' %in% tmp) {
        tmp[(length(tmp)-1):length(tmp)]  
      } else if ('FAILED' %in% tmp | 'RUNNING' %in% tmp) {
        c(NA,NA)
      }
    },USE.NAMES = F)
    if (is.list(a)){
      a = a[!is.na(a)]
      a = do.call(cbind,a)
    }
    a = a[,order(as.numeric(sub('K','',a[2,])))]
    a[1,] <- sapply(a[1,],timefunc) ## minutes
    a[2,] <- as.numeric(sub('K','',a[2,]))/1e6 ## GB
    
    a = data.frame(matrix(as.numeric(a),nrow=2), matrix(NA,2,4-ncol(a)), stringsAsFactors = FALSE)
    
    b <- data.frame(time = t(a[1,]), memory = t(a[2,]), numSamples = f, numCells = c(10,1e2,1e3,1e4), method = mtd, stringsAsFactors = FALSE)
    colnames(b) <- c('time', 'memory', 'numSamples', 'numCells', 'method')
    b
  })
  do.call(rbind,eff)
})
res <- do.call(rbind, res)
saveRDS(res,'/home-4/whou10@jhu.edu/work-zfs/whou10/metpred/scalability/result/time_memory.rds')
