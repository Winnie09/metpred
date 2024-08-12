rm(list=ls())
pdir <- '/home/whou10/data/whou/metpred/evaluate/eff/plot/'
source('/home/whou10/scratch16/whou10/resource/ggplot_theme.R')
library(ggplot2)
library(gridExtra)
library(reshape2)

## time
tb = readRDS('/home/whou10/data/whou/metpred/evaluate/eff/res/time_include_allcpg.rds')
tb[,1] = c('Chr10 CpGs', 'all CpGs')
t = matrix(as.numeric(tb[,-1]),nrow=nrow(tb))
dimnames(t) = list(tb[,1],colnames(tb)[-1])
t2 = melt(t)
colnames(t2) = c('type', 'numberSample', 'time')

## memory
m = readRDS('/home/whou10/data/whou/metpred/evaluate/eff/res/memory_include_allcpg.rds')
m[,1] = c('Chr10 CpGs', 'all CpGs')
rownames(m) = m[,1]
m = m[,-1]
dn = dimnames(m)
m = matrix(as.numeric(sub('K','',m)),ncol=ncol(m))
dimnames(m) = dn
m = m/(1024^2)
m2 = melt(m)
colnames(m2) = c('type', 'numberSample', 'memory')

## plot
p1 <- ggplot(data = t2, aes(x = numberSample, y = time, group = type,color = type)) +
  geom_point(size = 0.6) +
  geom_line(lwd = 0.2) +
  xlab(paste0("Number of samples")) +
  ylab('Time (min)') +
  theme(plot.margin = margin(5,5,5,5))
p2 <- ggplot(data = m2, aes(x = numberSample, y = memory, group = type, color = type)) +
  geom_point(size = 0.6) +
  geom_line(lwd = 0.2, color = 'grey') +
  xlab(paste0("Number of samples")) +
  ylab('Memory (GB)') +
  theme(plot.margin = margin(5,5,5,5))

pdf(paste0(pdir, 'eff_simu_wgbs_allcpg.pdf'), width = 5.5, height = 1.5)
grid.arrange(p1, p2, nrow = 1, padding = unit(0.01, "line") )
dev.off()

