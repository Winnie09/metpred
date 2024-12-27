## check tcga three DNA methylation tech: batch effect
## wgbs, 450k, epic
setwd('/home/whou10/data/whou/metpred/data/tcga/hg38/2024/combine/')
rdir = '/home/whou10/data/whou/metpred/data/tcga/hg38/2024/quality/res/'

## read in data
wgbs <- readRDS('wgbs.rds')  
array = readRDS('450k.rds')
epic = readRDS('EPIC.rds')

## use ramp: train on 450k predict on epic
## train 450k, predict on wgbs
int = intersect(rownames(wgbs), rownames(epic))
int = intersect(int, rownames(array))

set.seed(1)
cpgsel = sample(int, 1e3)
saveRDS(cpgsel, paste0(rdir, 'cpgsel1k.rds'))

