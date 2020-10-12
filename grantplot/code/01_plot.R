pddir <- '/home-4/whou10@jhu.edu/work-zfs/whou10/metpred/simu/plot/pd/'
pdir <- '/home-4/whou10@jhu.edu/work-zfs/whou10/metpred/grantplot/plot/'
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
mtdcolor <-  brewer.pal(8,'Set1')[-6]
names(mtdcolor) <- mtdorder


pdf(paste0(pdir, 'nullsimu.pdf'), width = 3.5, height = 2.5)
ggplot(data = pd, aes(x = Data, y = Num.FP, fill = Method, color = Data)) +
  geom_bar(stat = 'identity', position=position_dodge(), alpha = 0.5) +
  theme_classic() +
  scale_color_brewer(palette = 'Dark2') +
  # scale_fill_brewer(palette = 'Pastel1') +
  scale_fill_manual(values = mtdcolor) +
  ylab(bquote('Num.FP (' ~ log[10] ~ '-scaled )')) +
  coord_flip() + 
  theme(axis.text = element_text(size = 8, color = 'black')) +
  guides(fill= guide_legend("Method"), color = FALSE)
dev.off()


pd2 <- readRDS(paste0(pddir, 'perf_auc.rds')) ######### caution
pd2$method <- factor(as.character(pd2$method), levels = mtdorder)
pd2 <- pd2[!grepl('0.2.rds',pd2$type ), ]
pd2 <- pd2[!grepl('^10_',pd2$type ), ]

pdf(paste0(pdir, 'simu_auc.pdf'), width = 2.5, height = 2.6)
ggplot(data = pd2) + 
      geom_violin(aes(x = method, y = performance, fill = method), alpha = 0.2, scale = 'width')+
      geom_jitter(aes(x = method, y = performance, color = method),size = 0.2, alpha = 0.4, width = 0.2) +
      theme_classic()+
      theme(legend.position = 'none') +
      ylab('AUC')+ xlab('')+
      scale_fill_manual(values = mtdcolor)+
      scale_color_manual(values = mtdcolor)+
      coord_flip()+
      theme(axis.text = element_text(size = 8, color = 'black'))
dev.off()

ddir <- '/home-4/whou10@jhu.edu/work-zfs/whou10/metpred/simu/perf/sensfdr/200_200_0.2_0.8/'
sensfdr <- readRDS(paste0(ddir, 'limma_saver.rds.rds'))
pd3 <- data.frame(Real_FDR = sensfdr[,2], Reported_FDR = sensfdr[,3], stringsAsFactors = FALSE)
pdf(paste0(pdir, 'limma_fdr_curve.pdf'), width = 2.2, height = 2.2)
ggplot(data = pd3, aes(x = Reported_FDR, y = Real_FDR)) +
      geom_point(size = 0.5) +
      theme_classic()  +
      geom_abline(slope = 1, intercept = 0, color = 'red')+
      geom_vline(xintercept = 0.25, color = 'blue')+
      theme(axis.text = element_text(size = 8, color = 'black'))
dev.off()


ddir <- '/home-4/whou10@jhu.edu/work-zfs/whou10/metpred/scalability/result/'
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
time[,1] <- factor(time[,1], levels = rev(mtdorder))
time[,3] <- paste0(time[,3]*2, ' samples')
time[,3] <- factor(time[,3], levels = paste0(c(6,8,10,20), ' samples'))
time <- time[complete.cases(time),]
library(ggplot2)
library(RColorBrewer)

pdf(paste0(pdir, 'scalability.pdf'), width = 6, height = 2.4)
ggplot(time,aes(x=cn,y=time,col=met,group=met)) + 
  geom_point() +
  geom_line(size = 0.5) + 
  scale_color_manual(values = mtdcolor) +
  theme_classic() +
  scale_x_continuous(trans='log10') +
  ylab('Elapsed time (seconds)') + 
  xlab('Number of cells per sample')+
  ylim(c(0, 900)) +
  facet_wrap(~sn, nrow = 1) 
dev.off()


