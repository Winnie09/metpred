###########################################
## check tcga gene expression batch effect
## across those corresponding to three DNA methylation tech
setwd('/home/whou10/data/whou/metpred/data/tcga/hg38/2024/combine/')
rdir = '/home/whou10/data/whou/metpred/data/tcga/hg38/2024/quality/res/'
pdir = '/home/whou10/data/whou/metpred/data/tcga/hg38/2024/quality/plot/'

source('/home/whou10/data/whou/resource/myfunc/01_function.R')
ge = readRDS('ge.rds')
a <- PCA(ge, save.pca = T, plot.statistics=T, plot.dir = pdir, result.dir = rdir)
saveRDS(a, paste0(rdir, 'ge_pc.rds'))
