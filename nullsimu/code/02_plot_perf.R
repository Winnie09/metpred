perf <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/metpred/nullsimu/res/perf/num.fp.rds')
perf <- perf[!rownames(perf) %in% 'limmacell_saver', ]
perf <- log10(perf + 1)
rownames(perf) <- c('GLMM', 'limma', 'MAST', 'ourmethod', 'scDD', 't-test','wilcox')
library(reshape2)
library(ggplot2)
library(RColorBrewer)
pd = melt(perf)
colnames(pd) = c('Method','Data','Num.FP')
mtdorder = names(sort(tapply(pd[,'Num.FP'], list(pd[,'Method']), mean,na.rm=T), decreasing = T))
stat = tapply(pd[,'Num.FP'], list(pd[,'Method']), mean,na.rm=T)
pd$Method = factor(as.character(pd$Method), levels=mtdorder)
theme_hm <- function(method_vec) {
  theme_minimal() +
    theme(axis.text.x = element_text(size=12,color='black',angle=45,hjust=1),
          axis.text.y = element_text(size=12,color=ifelse(levels(method_vec)=='ourmethod','red','black')),
          plot.title=element_text(size=12),
          legend.title=element_text(size=12),
          legend.position = 'bottom')
}

pdir <- '/home-4/whou10@jhu.edu/scratch/Wenpin/metpred/nullsimu/plot/'
pdf(paste0(pdir, 'num.fp.pdf'), height = 5, width = 3.5)
ggplot() + geom_tile(data=pd,aes(x=Data,y=Method,fill=Num.FP)) + theme_hm(pd$Method) +
  scale_fill_gradientn(colors=rev(brewer.pal(9,'YlGnBu')))  + 
  xlab('') + ylab('')
dev.off()
