library(ggplot2)
library(reshape2)
res <- sapply(c(1,5,10,20),function(i) {
  readRDS(paste0('/home-4/zji4@jhu.edu/scratch/metpred/cellline/predict/singlesite/sc/cor_',i,'.rds'))[[1]]
})

d <- readRDS(paste0('/home-4/zji4@jhu.edu/scratch/metpred/cellline/predict/singlesite/bulk/cor.rds'))
res <- cbind(res,d[[1]])
row.names(res) <- c('A549','GM12878','IMR-90','K562')
colnames(res) <- c(paste0('sc_',c(1,5,10,20)),'bulk')
pd <- melt(res)
td <- apply(res,2,median,na.rm=T)
td <- data.frame(method=names(td),correlation=round(td,3),stringsAsFactors = F)
td$pos <- td$correlation + 0.1
colnames(pd) <- c('celltype','method','correlation')
pdf(paste0('/home-4/zji4@jhu.edu/scratch/metpred/cellline/plot/acrosssite.pdf'))
print(ggplot() + geom_boxplot(data=pd,aes(x=method,y=correlation)) + geom_point(data=pd,aes(x=method,y=correlation,col=celltype)) + geom_text(data=td,aes(x=method,y=pos,label=correlation)) + theme_classic())
dev.off()  


