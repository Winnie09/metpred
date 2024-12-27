## check tcga three DNA methylation tech: batch effect
## wgbs, 450k, epic
setwd('/home/whou10/data/whou/metpred/data/tcga/hg38/2024/combine/')
rdir = '/home/whou10/data/whou/metpred/data/tcga/hg38/2024/quality/res/'

library(matrixStats)

## read in data
array = readRDS('450k.rds')
epic = readRDS('EPIC.rds')
ge = readRDS('ge.rds')

## use ramp: train on 450k predict on epic
## train 450k, predict on wgbs
cpgsel = readRDS(paste0(rdir, 'cpgsel1k.rds'))

trainexpr = ge[, colnames(array)]
trainexpr = trainexpr[rowMeans(ge) > 1, ]
trainexpr = trainexpr[rowSds(ge) > 0.2, ]
print(str(trainexpr))


source('/home/whou10/data/whou/metpred/software/trainpredict.R')
predepic <- trainpredict(trainexpr = trainexpr,
             trainmeth = array[cpgsel, ],
             testexpr = ge[, colnames(epic)],
             clunumlist = 1e3, lambdalist = 0.1)
saveRDS(predepic, paste0(rdir, 'trainarray_predepic.rds'))

