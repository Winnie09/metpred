## calculate pbmc perf
rm(list=ls())
method = commandArgs(trailingOnly = T)[1] ##. neargene or permu
# key = 'pred_var_10000_CpG'
# key = 'pred_varmedian_0673446_CpG'
key = 'allcpg'
print(method)

source('/home/whou10/scratch16/whou10/resource/myfunc/01_function.R')
ddir <- '/home/whou10/data/whou/metpred/data/pbmc/imp/res/'
rdir = '/home/whou10/data/whou/metpred/evaluate/pbmc/predsc/imp/'

## load prediction
pred = readRDS(paste0(rdir, '/competing/combine/', method, '.rds'))

## load gs, match dim with pred
gs = readRDS('/home/whou10/data/whou/metpred/evaluate/pbmc/predsc/imp/res/eval/eval_goldstandard_ctmean.rds')
o = intersect(rownames(pred), rownames(gs))
pred = pred[o, ]
gs = gs[o, ]
pred = pred[, colnames(gs)]

## across-sample PCC
pd = acrossRowCor_plotdata(pred = pred, goldstandard = gs)
saveRDS(pd, paste0(rdir, 'perf/', method, '/', key, '_acrosssample_pcc.rds'))

## across-cpg PCC
id = complete.cases(pred) & complete.cases(gs)
pd2 = acrossRowCor_plotdata(pred = t(pred[id,]), 
                            goldstandard = t(gs[id,]))
saveRDS(pd2, paste0(rdir, 'perf/', method, '/', key, '_acrosscpg_pcc.rds'))
