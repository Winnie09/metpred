res <- readRDS('/home-4/whou10@jhu.edu/work-zfs/whou10/metpred/scalability/result/time_memory.rds')
pdir <- '/home-4/whou10@jhu.edu/work-zfs/whou10/metpred/scalability/plot/'
res$numSamples <- factor(res$numSamples, levels = c('5', '10'))
library(ggplot2)
library(RColorBrewer)
library(reshape2)
pdf(paste0(pdir, 'time.pdf'), width = 6, height = 2.4)
ggplot(res, aes(x=numCells,y=time,col=method,group=method)) + 
  geom_point() +
  geom_line(size = 0.5) + 
  scale_color_brewer(palette = 'Set1') +
  theme_classic() +
  scale_x_continuous(trans='log10') +
  ylab('Elapsed time (minutes)') + 
  xlab('Number of cells per sample')+
  facet_wrap(~numSamples, nrow = 1) 
dev.off()

pdf(paste0(pdir, 'memory.pdf'), width = 6, height = 2.4)
ggplot(res, aes(x=numCells,y=memory,col=method,group=method)) + 
  geom_point() +
  geom_line(size = 0.5) + 
  scale_color_brewer(palette = 'Set1') +
  theme_classic() +
  scale_x_continuous(trans='log10') +
  ylab('Memory (MaxRSS, GB)') + 
  xlab('Number of cells per sample')+
  facet_wrap(~numSamples, nrow = 1) 
dev.off()

