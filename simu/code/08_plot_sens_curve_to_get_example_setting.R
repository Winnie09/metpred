ddir <- '/home-4/whou10@jhu.edu/work-zfs/whou10/metpred/simu/perf/sensfdr/'
## 200_200_0.2_0.8
af <- list.files(ddir)
af <- af[af != '2_2_0.4_0.4']
pdf("/home-4/whou10@jhu.edu/work-zfs/whou10/metpred/simu/plot/plot/sensfdr_curve_limma.pdf", width = 30, height = 24)
par(mfrow = c(15,12))
for (f in af) {
  print(f)
  sensfdr <- readRDS(paste0(ddir, f, '/limma_saver.rds.rds'))
  plot(sensfdr[,1] ~ sensfdr[,2], ylab='Sensitivity', xlab='RealFDR',pch=20)
  abline(a=0,b=1,col='red')
  plot(sensfdr[,2] ~ sensfdr[,3], ylab='RealFDR', xlab='ReportedFDR',pch=20)
  abline(a=0,b=1,col='red')
}
dev.off()


library(ggplot2)
library(RColorBrewer)
library(gridExtra)
pdf('/home-4/whou10@jhu.edu/work-zfs/whou10/metpred/simu/plot/sens_fdr_curve.pdf', width = 31.5, height = 40)
plist <- list()
afd <- setdiff(list.files('/home-4/whou10@jhu.edu/work-zfs/whou10/metpred/simu/perf/sensfdr/'), c('2_2_0.4_0.4'))
for (fd in afd){
  print(fd)
  ddir <- paste0('/home-4/whou10@jhu.edu/work-zfs/whou10/metpred/simu/perf/sensfdr/', fd, '/')## 200_200_0.2_0.8
  af <- list.files(ddir)
  allsens <- lapply(af, function(f){
    m <- gsub('_.*','',sub('.rds.rds','',f))
    if (m == 'glmm') m <- 'GLMM'
    if (m == 'wilcoxon') m <- 'wilcox'
    if (m == 't') m <- 't-test'
    data.frame(readRDS(paste0(ddir, f)), Method = m)
  })
  pd3 <- do.call(rbind, allsens)
  pd3[,4] <- factor(pd3[,4])
  plist[[fd]] <- ggplot(data = pd3, aes(x = Real_FDR, y = Sensitivity, color=Method, group = Method)) +
    geom_line(size = 1) +
    theme_classic()  +
    geom_abline(slope = 1, intercept = 0, color = 'red')+
    # geom_vline(xintercept = 0.25, color = 'blue')+
    theme(axis.text = element_text(size = 8, color = 'black')) +
    xlim(c(0,0.1)) + 
    ggtitle(fd) + 
    scale_color_brewer(palette='Set1')
}
grid.arrange(grobs = plist, nrow = 16)  
dev.off()



library(ggplot2)
library(RColorBrewer)
library(gridExtra)
pdf('/home-4/whou10@jhu.edu/work-zfs/whou10/metpred/simu/plot/sens_fdr_curve_selected.pdf', width = 13, height = 7.5)
plist <- list()
afd <- c('2_10_0.2_0.3',
'2_10_0.4_0.3',
'2_2_0.2_0.3',
'2_2_0.4_0.3',
'4_200_0.3_0.4',
'4_2_0.4_0.4',
'4_2_0.2_0.4',
'4_100_0.3_0.4',
'2_2_0.3_0.3',
'2_10_0.3_0.3')
for (fd in afd){
  print(fd)
  ddir <- paste0('/home-4/whou10@jhu.edu/work-zfs/whou10/metpred/simu/perf/sensfdr/', fd, '/')
  af <- list.files(ddir)
  allsens <- lapply(af, function(f){
    m <- gsub('_.*','',sub('.rds.rds','',f))
    if (m == 'glmm') m <- 'GLMM'
    if (m == 'wilcoxon') m <- 'wilcox'
    if (m == 't') m <- 't-test'
    data.frame(readRDS(paste0(ddir, f)), Method = m)
  })
  pd3 <- do.call(rbind, allsens)
  pd3[,4] <- factor(pd3[,4])
  plist[[fd]] <- ggplot(data = pd3, aes(x = Real_FDR, y = Sensitivity, color=Method, group = Method)) +
    geom_line(size = 1) +
    theme_classic()  +
    geom_abline(slope = 1, intercept = 0, color = 'red')+
    # geom_vline(xintercept = 0.25, color = 'blue')+
    theme(axis.text = element_text(size = 8, color = 'black')) +
    xlim(c(0,0.05)) + 
    ggtitle(fd) + 
    scale_color_brewer(palette='Set1')
}
grid.arrange(grobs = plist, nrow = 4)  
dev.off()


# ----------------
library(ggplot2)
library(RColorBrewer)
library(gridExtra)
# pdf('/home-4/whou10@jhu.edu/work-zfs/whou10/metpred/simu/plot/sens_fdr_curve.pdf', width = 31.5, height = 40)
plist <- list()
afd <- setdiff(list.files('/home-4/whou10@jhu.edu/work-zfs/whou10/metpred/simu/perf/sensfdr/'), c('2_2_0.4_0.4'))
## sensitivity at real fdr 0.05
mat <- sapply(afd, function(fd){
  print(fd)
  ddir <- paste0('/home-4/whou10@jhu.edu/work-zfs/whou10/metpred/simu/perf/sensfdr/', fd, '/')## 200_200_0.2_0.8
  af <- list.files(ddir)
  allsens <- sapply(af, function(f){
    m <- gsub('_.*','',sub('.rds.rds','',f))
    if (m == 'glmm') m <- 'GLMM'
    if (m == 'wilcoxon') m <- 'wilcox'
    if (m == 't') m <- 't-test'
    tmp <- data.frame(readRDS(paste0(ddir, f)), Method = m)
    tmp[which(tmp[,2]>0.05)[1], 1]
  })
  names(allsens) <- gsub('.rds.rds', '', sub('_.*', '', names(allsens)))
  allsens
})  

library(pheatmap)
pheatmap(t(mat[-c(1,5,6,7), ]), cluster_cols = F, cluster_rows = F,
         color = rainbow(100))


