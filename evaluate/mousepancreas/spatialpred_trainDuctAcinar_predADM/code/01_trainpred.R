s =  as.character(commandArgs(trailingOnly = T)[1]) # s = 'A46_F_P_K4_D4'; s = 'A47_F_P_K4_D7'
source('/home/whou10/scratch16/whou10/resource/startup.R')
setwd('/home/whou10/data/whou/metpred/')
dir.r <- paste0('evaluate/mousepancreas/spatialpred_trainDuctAcinar_predADM/res/', s)
dir.create(dir.r, recursive = T, showWarnings = F)
cpglist <- readRDS(paste0('evaluate/mousepancreas/celltypespecific_cv/cpg_group/', s, '.rds'))
id = as.numeric(commandArgs(trailingOnly = T)[2]) ## id = 1
print(id)
id = length(cpglist) - id + 1 ## run largest var CpG group first
if (id <=0 ) {stop('cpglist length smaller than $j.')}
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
intct = sub('.*:','',int)
intsamp <- sub(':.*','',int)
int2 = int[intct %in% c('acinar', 'duct')]
samp <- sub(':.*','',int2)
trainid <- int2[which(samp != s)]
testid <- setdiff(int[intsamp == s], trainid)

library(Matrix)
expr <- expr[, sub(':.*','',colnames(expr)) == s]
source('/home/whou10/data/whou/metpred/software/trainpredict.R')
source('/home/whou10/scratch16/whou10/resource/startup.R')
pred <- trainpredict(trainexpr=r[,trainid],
                     testexpr=expr,
                     trainmeth=w[cpglist[[id]], trainid])
saveRDS(pred, file_name)



