setwd('/home-4/whou10@jhu.edu/work-zfs/whou10/metpred/scalability/')
time <- readRDS('./result/time.rds')
library(ggplot2)
library(RColorBrewer)
pdf('./plot/time.pdf', width = 4, height = 2.8)
ggplot(time,aes(x=cn,y=time,col=met,group=met)) + 
  geom_point(alpha = 0.5) +
  geom_line(size = 0.5) + 
  scale_color_brewer(palette = 'Set2') +
  theme_classic() +
  scale_x_continuous(trans='log10') +
  ylab('Elapsed time (seconds)') + 
  xlab('Number of cells per sample')
dev.off()
