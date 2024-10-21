rm(list=ls())

rdir = '/home/whou10/data/whou/metpred/evaluate/pbmc/predsc/imp/perf/ramp/'
# d = readRDS('/home/whou10/data/whou/rna_imputation/result/procimpute/pbmc/saver/sorted.rds')
d = readRDS('/home/whou10/data/whou/metpred/data/pbmc/imp/res/rna_test_sc_sub300.rds')
str(d)
source('/home/whou10/scratch16/whou10/resource/myfunc/01_function.R')
pr = pca_lmFilter(genebycellmat = d, topnum = 2e3)
saveRDS(pr, paste0(rdir, 'pr.rds'))

u = UMAP(pr$x)
saveRDS(u, paste0(rdir, 'umap.rds'))

## ====================
## plot
## =====================
pdir = '/home/whou10/data/whou/metpred/evaluate/pbmc/predsc/imp/plot/'
u = readRDS(paste0(rdir, 'umap.rds'))
pd = data.frame(u1 = u[,1], u2 = u[,2], 
                celltype = sub(':.*', '', rownames(u)))
library(ggplot2)
library(RColorBrewer)

pdf(paste0(pdir, 'umap_sc_sub300.pdf'), width = 1.8, height = 1.8)
plotumap(plotdata = pd)
dev.off()
