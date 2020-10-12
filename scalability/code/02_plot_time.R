ddir <- '/home-4/whou10@jhu.edu/work-zfs/whou10/metpred/scalability/result/'
pdir <- '/home-4/whou10@jhu.edu/work-zfs/whou10/metpred/scalability/plot/'
af <- list.files(ddir)
f = af[1]
tmp <- lapply(af, function(f){
  return(readRDS(paste0(ddir, f)))
})
time <- do.call(rbind, tmp)  

time[,1] <- sub('_saver','',time[,1])
time[which(time[,1] == 'wilcoxon'),1] <- 'wilcox'
time[which(time[,1] == 'glmm'),1] <- 'GLMM'
time[which(time[,1] == 't'),1] <- 't-test'
time[which(time[,1] == 'RAISIN'),1] <- 'ourmethod'
time[,3] <- paste0(time[,3]*2, ' samples')
time[,3] <- factor(time[,3], levels = paste0(c(6,8,10,20), ' samples'))
time <- time[complete.cases(time),]
library(ggplot2)
library(RColorBrewer)

pdf(paste0(pdir, 'scalability.pdf'), width = 6, height = 2.4)
ggplot(time,aes(x=cn,y=time,col=met,group=met)) + 
  geom_point() +
  geom_line(size = 0.5) + 
  scale_color_brewer(palette = 'Set1') +
  theme_classic() +
  scale_x_continuous(trans='log10') +
  ylab('Elapsed time (seconds)') + 
  xlab('Number of cells per sample')+
  facet_wrap(~sn, nrow = 1) 
dev.off()

