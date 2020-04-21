library(ggplot2)
library(reshape2)
for (sdc in c(0,0.1,0.2)) {
  res <- sapply(c(1,5,10,20),function(i) {
    d <- readRDS(paste0('/home-4/zji4@jhu.edu/scratch/metpred/cellline/predict/singlesite/sc/cor_',i,'.rds'))
    d[[2]][d[[3]] > sdc]
  })
  
  d <- readRDS(paste0('/home-4/zji4@jhu.edu/scratch/metpred/cellline/predict/singlesite/bulk/cor.rds'))
  res <- cbind(res,d[[2]][d[[3]] > sdc])
  colnames(res) <- c(paste0('sc_',c(1,5,10,20)),'bulk')
  pd <- melt(res)
  td <- apply(res,2,median,na.rm=T)
  td <- data.frame(method=names(td),correlation=round(td,3),stringsAsFactors = F)
  td$pos <- td$correlation + 0.1
  colnames(pd) <- c('site','method','correlation')
  pdf(paste0('/home-4/zji4@jhu.edu/scratch/metpred/cellline/plot/acrosssample_sd',sdc,'.pdf'))
  print(ggplot() + geom_boxplot(data=pd,aes(x=method,y=correlation)) + geom_text(data=td,aes(x=method,y=pos,label=correlation)) + theme_classic())
  dev.off()  
}

