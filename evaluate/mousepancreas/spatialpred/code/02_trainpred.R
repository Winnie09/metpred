## sd of 
setwd('/home/whou10/data/whou/metpred/evaluate/mousepancreas/')
dir.r <- 'spatialpred/res/'
source('/home/whou10/scratch16/whou10/resource/startup.R')

expr <- readRDS('/home/whou10/data/whou/metpred/data/mousepancreas/spatial/normexpr.rds')
alls <- sub(':.*', '', colnames(expr))
sid <- as.numeric(commandArgs(trailingOnly = T)[[1]])
print(sid)
s = unique(alls)[sid]

## load cell type specific expression and wgbs data
r <- readRDS('/home/whou10/data/whou/metpred/data/mousepancreas/spatial/cse.rds')
w <- readRDS('/home/whou10/data/whou/metpred/data/mousepancreas/wgbs/bs.rds')
genecpg <- readRDS(paste0('spatialpred/genecpg/', s, '.rds'))
gene.sel <- unique(genecpg[,2])

## train model and predict
int <- intersect(colnames(r),colnames(w))
samp <- sub(':.*','',int)
trainid <- int[which(samp != s)]
testid <- setdiff(int, trainid)
library(Matrix)
se <- sub(':.*','',colnames(expr))
expr <- expr[, se == s]
str(expr)
source('/home/whou10/data/whou/metpred/software/trainpredict.R')
source('/home/whou10/scratch16/whou10/resource/startup.R')
pred <- trainpredict(trainexpr=r[!rownames(r) %in% gene.sel,trainid],
                     testexpr=expr[!rownames(expr) %in% gene.sel,],
                     trainmeth=w[genecpg[,1], trainid])
saveRDS(pred, paste0('spatialpred/res/pred_', s, '.rds'))

