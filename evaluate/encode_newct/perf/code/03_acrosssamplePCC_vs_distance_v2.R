rm(list=ls())
## acrosssample pcc vs distance
## each dot is a cv

## load pcc
library(ggplot2)
pal <- c('Ramp'='red','Nearest gene'='orange','Permute'='royalblue')
pd <- NULL

cvid = 1
pdlist <- mpdlist <- list()
for (cvid in 1:10){
  type = 'acrosssample'
  rdir <- '/home/whou10/data/whou/metpred/evaluate/encode_newct/perf/res/'
  p1 <- data.frame(method='Nearest gene',type=type,pcc=median(readRDS(file=paste0(rdir, 'neargene/', type, '/cv', cvid, '.rds'))$cor, na.rm = T),stringsAsFactors = F)
  p2 <- data.frame(method='Permute',type=type,pcc=median(readRDS(paste0(rdir, 'permu/', type, '/cv', cvid, '.rds'))$cor,na.rm=T),stringsAsFactors = F)
  p3 <- data.frame(method='Ramp',type=type,pcc=median(readRDS(paste0(rdir, 'ramp/', type, '/cv', cvid, '.rds'))$cor,na.rm=T),stringsAsFactors = F)
  pdlist[[cvid]] = data.frame(rbind(p1,p2,p3), cv=cvid)
}

pd = do.call(rbind, pdlist)

## load rna-seq distance
ddir = '/home/whou10/data/whou/metpred/evaluate/encode_newct/pc/rna_cosineSimilarity/'
cvid = 1
disv <- sapply(1:10, function(cvid){
  dis = 1-readRDS(paste0(ddir, cvid, '.rds'))
  median(apply(dis, 1, min, na.rm=T), na.rm=T) ## median of minimum
  
})
names(disv) = 1:10

pd2 = cbind(pd, distance=disv[pd$cv])

library(ggplot2)

pal <- c('Ramp'='red','Nearest gene'='orange','Permute'='royalblue')
pdir = '/home/whou10/data/whou/metpred/evaluate/encode_newct/perf/plot/compare/'
saveRDS(pd2, paste0(pdir, 'acrosssamplePCC_vs_distance_v2.rds'))

## each data point is, in a cross validation, the median of minimum of the testing cell type to training cell types 
pdf(paste0(pdir, '/acrosssamplePCC_vs_distance_v2.pdf'), width=3.2,height=2.2)
ggplot(pd2[pd2$method=='Ramp',],aes(x=distance,y=pcc)) + 
  geom_point(size = 1, stroke = 0, width = 0.2) + 
  geom_smooth(method = "lm", se = F) +
  theme_classic() + 
  theme(legend.position = 'bottom') + 
  xlab('Distance between training and testing') + 
  ylab('Across-sample PCC') + 
#  scale_color_manual(values=pal) + 
  theme(axis.text.x = element_text(angle = 30,vjust=0.5),legend.title = element_blank()) 
# scale_y_continuous(breaks = c(-0.25,0, 0.25, 0.5,0.75))
dev.off()


