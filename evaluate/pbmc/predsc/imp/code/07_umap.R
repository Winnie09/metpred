rm(list=ls())
d = readRDS('/home/whou10/data/whou/rna_imputation/result/procimpute/pbmc/saver/sorted.rds')
str(d)
source('/home/whou10/scratch16/whou10/resource/myfunc/01_function.R')

pr = pca_lmFilter(genebycellmat = d, topnum = 2e3)
str(pr)
u = UMAP(pr)
str(u)
pdir = '/home/whou10/data/whou/metpred/evaluate/pbmc/predsc/imp/plot/'
rdir = '/home/whou10/data/whou/metpred/evaluate/pbmc/predsc/imp/res/'
saveRDS(u, paste0(rdir, 'umap.rds'))



  



