## check tcga three DNA methylation tech: batch effect
## wgbs, 450k, epic
setwd('/home/whou10/data/whou/metpred/data/tcga/hg38/2024/combine/')
rdir = '/home/whou10/data/whou/metpred/data/tcga/hg38/2024/quality/res/'
library(matrixStats)
## read in data
wgbs = readRDS('wgbs.rds')
array = readRDS('450k.rds')
ge = readRDS('ge.rds')

## use ramp: train on 450k predict on epic
## train 450k, predict on wgbs
cpgsel = readRDS(paste0(rdir, 'cpgsel1k.rds'))

trainexpr = ge[, colnames(wgbs)]
trainexpr = trainexpr[rowMeans(ge) > 1, ]
trainexpr = trainexpr[rowSds(ge) > 0.2, ]
print(str(trainexpr))


source('/home/whou10/data/whou/metpred/software/trainpredict.R')
predwgbs <- trainpredict(trainexpr = trainexpr,
             trainmeth = array[cpgsel, ],
             testexpr = ge[, colnames(wgbs)],
             clunumlist = 1e3, lambdalist = 0.1)
saveRDS(predwgbs, paste0(rdir, 'trainarray_predwgbs.rds'))


