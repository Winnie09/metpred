rdir = '/home/whou10/data/whou/metpred/evaluate/pbmc/predsc/imp/perf/ramp/'
source('/home/whou10/scratch16/whou10/resource/myfunc/01_function.R')
d = readRDS(paste0(rdir, 'predicted_cpg_sd0.3.rds'))
pr = pca_lmFilter(genebycellmat = d, topnum = 2e3)
saveRDS(pr, paste0(rdir, 'pred_dnam_pr.rds'))

u = UMAP(pr$x)
saveRDS(u, paste0(rdir, 'pred_dnam_umap.rds'))

## ====================
## plot
## =====================
pdir = '/home/whou10/data/whou/metpred/evaluate/pbmc/predsc/imp/plot/'
u = readRDS(paste0(rdir, 'pred_dnam_umap.rds'))
pd = data.frame(u1 = u[,1], u2 = u[,2], 
                celltype = sub(':.*', '', rownames(u)))
# pd[which(pd[,3]=='cd14_monocytes'),3] = 'monocyte'
library(ggplot2)
library(RColorBrewer)
pdf(paste0(pdir, 'umap_pred_dnam.pdf'), width = 1.8, height = 1.8)
plotumap(plotdata = pd)
dev.off()

rm(list=ls())
