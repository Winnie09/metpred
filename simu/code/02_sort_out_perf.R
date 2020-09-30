rm(list=ls())
source('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/function/01_function.R')
ddir1 <- '/home-4/whou10@jhu.edu/scratch/Wenpin/metpred/simu/res/'
ddir2 <- '/home-4/whou10@jhu.edu/scratch/Wenpin/metpred/simu/data/data/'
rdir1 <- '/home-4/whou10@jhu.edu/scratch/Wenpin/metpred/simu/perf/sensfdr/'
rdir2 <- '/home-4/whou10@jhu.edu/scratch/Wenpin/metpred/simu/perf/perf/'
pdir <- '/home-4/whou10@jhu.edu/scratch/Wenpin/metpred/simu/plot/'
af <- list.files(ddir1)

library(parallel)
tmp <- mclapply(af, function(f){
  print(f)
  dir.create(paste0(rdir1, f))
  diffgn <- readRDS(paste0(ddir2, f, '/diffgn.rds'))
  am <- list.files(paste0(ddir1, f))
  perf <- t(sapply(am, function(m){
    res <- readRDS(paste0(ddir1, f, '/', m))
    sensfdr <- SensFdr(diffgn, res)
    saveRDS(sensfdr, paste0(rdir1, f, '/', m, '.rds'))
    AreaUnderSensFdr(sensfdr)
  }))
  saveRDS(perf, paste0(rdir2, f, '.rds'))
  return(0)
}, mc.cores = detectCores())

# -----------
# plot result
# ------------
af <- list.files(rdir2)
af <- af[!grepl('0.2.rds', af)]
perf <- t(sapply(af, function(f){
  tmp <- readRDS(paste0(rdir2, f))
  rownames(tmp) <- sub('.rds','',rownames(tmp))
  tmp[,1]
}))
colnames(perf) <- c('GLMM', 'limma', 'MAST', 'ourmethod', 'scDD', 't-test','wilcox')
library(reshape2)
pd1 <- melt(perf)
colnames(pd1) <- c('type', 'method', 'performance')
pd1 <- cbind(pd1, probRead = gsub('.*_','',sub('.rds','',pd1$type)))

perf <- t(sapply(af, function(f){
  tmp <- readRDS(paste0(rdir2, f))
  rownames(tmp) <- sub('.rds','',rownames(tmp))
  tmp[,2]
}))
colnames(perf) <- c('GLMM', 'limma', 'MAST', 'ourmethod', 'scDD', 't-test','wilcox')
library(reshape2)
pd2 <- melt(perf)
colnames(pd2) <- c('type', 'method', 'performance')
pd2 <- cbind(pd2, probRead = gsub('.*_','',sub('.rds','',pd2$type)))
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

#  check detailed performance w.r.t signal strength
p1 <- ggplot(data = pd1, aes(x = method, y = performance, fill = probRead)) + 
  geom_boxplot(alpha = 0.2)+
  theme_classic() +
  theme(legend.position = 'none') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylab('FDR.Diff')
p2 <- ggplot(data = pd2, aes(x = method, y = performance, color = probRead)) + 
  geom_boxplot()+
  theme_classic()+
  theme(legend.position = 'none') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylab('AUC')
pdf(paste0(pdir, 'perf_signal.pdf'), width = 10, height = 4)
grid.arrange(p1, p2, nrow = 1)
dev.off()

