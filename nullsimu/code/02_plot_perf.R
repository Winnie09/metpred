perf <- readRDS('/home-4/whou10@jhu.edu/work-zfs/whou10/metpred/nullsimu/res/perf/num.fp.rds')
perf <- perf[!rownames(perf) %in% 'limmacell_saver', ]
perf <- log10(perf + 1)
rownames(perf) <- c('GLMM', 'limma', 'MAST', 'ourmethod', 'scDD', 't-test','wilcox')
library(reshape2)
library(ggplot2)
library(RColorBrewer)
pd = melt(perf)
colnames(pd) = c('Method','Data','Num.FP')
mtdorder = rev(names(sort(tapply(pd[,'Num.FP'], list(pd[,'Method']), mean,na.rm=T), decreasing = T)))
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

pdir <- '/home-4/whou10@jhu.edu/work-zfs/whou10/metpred/nullsimu/plot/'
pdf(paste0(pdir, 'num.fp.pdf'), height = 4.5, width = 3)
ggplot() + geom_tile(data=pd,aes(x=Data,y=Method,fill=Num.FP)) + theme_hm(pd$Method) +
  scale_fill_gradientn(colors=rev(brewer.pal(9,'YlGnBu')))  + 
  xlab('') + ylab('')
dev.off()

pdf(paste0(pdir, 'num_fp_voilin.pdf'), height = 4, width = 6)
ggplot(data=pd) + 
  geom_violin(aes(x=Method, y=Num.FP, fill = Method), alpha = 0.2, scale = 'width') + 
  geom_jitter(aes(x = Method, y = Num.FP, color = Data), width=0.25, alpha=0.5) +
  xlab('') + ylab('log10(Num.FP+1)')+
  theme_classic()
dev.off()

pdf(paste0(pdir, 'num_fp_barplot.pdf'), height = 4, width = 4)
ggplot(data = pd, aes(x = Data, y = Num.FP, fill = Method, color = Data)) +
  geom_bar(stat = 'identity', position=position_dodge()) +
  theme_classic() +
  scale_color_brewer(palette = 'Dark2') +
  scale_fill_brewer(palette = 'Pastel1') +
  ylab(bquote('Num.FP (' ~ log[10] ~ '-scaled )')) +
  coord_flip() 
dev.off()
