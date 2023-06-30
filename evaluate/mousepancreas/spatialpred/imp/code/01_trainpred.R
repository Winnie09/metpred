 s =  as.character(commandArgs(trailingOnly = T)[[1]]) # s = 'A46_F_P_K4_D4'
source('/home/whou10/scratch16/whou10/resource/startup.R')
setwd('/home/whou10/data/whou/metpred/')
dir.create(dir.r <- paste0('evaluate/mousepancreas/spatialpred/imp/res/', s), recursive = T)
cpglist <- readRDS(paste0('evaluate/mousepancreas/celltypespecific_cv/cpg_group/', s, '.rds'))
id = as.numeric(commandArgs(trailingOnly = T)[[2]])
id = length(cpglist) - id + 1 ## run largest var CpG group first
file_name <- paste0(dir.r, '/pred_cpg_sd_', names(cpglist)[id], '.rds')
if (file.exists(file_name)) {
  stop("File already exists. Stopping execution.")
}

## load imputed visium expr, cell type specific expression and wgbs data
expr <- readRDS(file="data/mousepancreas/spatial/procimpute/res/all/imputednormexpr.rds")
r <- readRDS('data/mousepancreas/spatial/imputed_nnls/res/all/cse.rds')
w <- readRDS('data/mousepancreas/wgbs/bs.rds')

## train model and predict
int <- intersect(colnames(r),colnames(w))
samp <- sub(':.*','',int)
trainid <- int[which(samp != s)]
testid <- setdiff(int, trainid)
library(Matrix)
expr <- expr[, sub(':.*','',colnames(expr)) == s]
source('/home/whou10/data/whou/metpred/software/trainpredict.R')
source('/home/whou10/scratch16/whou10/resource/startup.R')
pred <- trainpredict(trainexpr=r[,trainid],
                     testexpr=expr,
                     trainmeth=w[cpglist[[id]], trainid])
saveRDS(pred, file_name)

