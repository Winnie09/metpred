rm(list=ls())
setwd('/home/whou10/data/whou/metpred/')
df = read.table('./evaluate/eff/res/efficiency_detail_scores.csv', sep = ',', header = T)
str(df)
df
pdir <- '/home/whou10/data/whou/metpred/evaluate/eff/plot/'
source('/home/whou10/scratch16/whou10/resource/ggplot_theme.R')

library(ggplot2)
library(gridExtra)
pdf(paste0(pdir, 'eff_simu_wgbs_chr10.pdf'), width = 5, height = 2.5)
p1 <- ggplot(data = df, aes(x = numberSample, y = t)) +
  geom_point(size = 0.6) +
  geom_line(lwd = 0.2, color = 'grey') +
  xlab(paste0("Number of samples")) +
  ylab('Time (min)')+
  ggtitle(paste0('Scalability: ', round(unique(df$scalability), 3)))


p2 <- ggplot(data = df, aes(x = numberSample, y = m)) +
  geom_point(size = 0.6) +
  geom_line(lwd = 0.2, color = 'grey') +
  xlab(paste0("Number of samples")) +
  ylab('Memory (GB)') 

grid.arrange(p1, p2, nrow = 1)
dev.off()
