pdir <- '/home-4/whou10@jhu.edu/work-zfs/whou10/metpred/simu/plot/plot/'
pddir <- '/home-4/whou10@jhu.edu/work-zfs/whou10/metpred/simu/plot/pd/'

pd1 <- readRDS(paste0(pddir, 'perf_fdrdiff.rds')) ######### caution
pd2 <- readRDS(paste0(pddir, 'perf_auc.rds')) ######### caution

library(ggplot2)
p1 <- ggplot(data = pd1, aes(x = method, y = performance, fill = method)) + 
  geom_violin(alpha = 0.3, scale = 'width')+
  theme_classic() +
  theme(legend.position = 'none') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylab('FDR.Diff') + xlab('')+
  scale_fill_brewer(palette = 'Set1')
p2 <- ggplot(data = pd2, aes(x = method, y = performance, fill = method)) + 
  geom_violin(alpha = 0.3, scale = 'width')+
  theme_classic()+
  theme(legend.position = 'none') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylab('AUC')+ xlab('')+
  scale_fill_brewer(palette = 'Set1')
library(gridExtra)
pdf(paste0(pdir, 'perf_method.pdf'), width = 5, height = 2.2)
grid.arrange(p1, p2, nrow = 1)
dev.off()
