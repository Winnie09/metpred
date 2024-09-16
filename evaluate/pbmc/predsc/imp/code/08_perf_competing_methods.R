## calculate pbmc perf
rm(list=ls())
method = commandArgs(trailingOnly = T)[1]
key = 'pred_var_10000_CpG'
print(method)

# key = 'pred_varmedian_0673446_CpG'

source('/home/whou10/scratch16/whou10/resource/myfunc/01_function.R')
ddir <- '/home/whou10/data/whou/metpred/data/pbmc/imp/res/'
w.te = readRDS(paste0(ddir, 'wgbs_test.rds'))

rdir = '/home/whou10/data/whou/metpred/evaluate/pbmc/predsc/imp/'
# pred = readRDS(paste0(rdir, 'res/pred/pred_varmedian_30673446_CpG.rds'))
# pred = readRDS(paste0(rdir, 'res/pred/pred_allcpg.rds'))
pred = readRDS(paste0(rdir, 'res/pred/pred_var_100000_CpG.rds'))

res = readRDS('/home/whou10/data/whou/metpred/evaluate/pbmc/predsc/imp/competing/combine/neargene.rds')
pred.agg = res[rownames(pred), ]

w.te2 = w.te[rownames(pred),]

# num [1:1644, 1:7] 0.816 0.769 0.89 0.298 0.92 ...
# - attr(*, "dimnames")=List of 2
# ..$ : chr [1:1644] "chr22_46144556" "chr10_30482068" "chr18_79422099" "chr2_223918843" ...
# ..$ : chr [1:7] "b_cells" "cd14_monocytes" "cd56_nk" "cytotoxic_t" ...
w.te.agg = aggregatefunc2(d = w.te2, 
                          by = c('cd56_nk', rep('cd14_monocytes', 3), 'b_cells', rep('t_cells', 4)), 
                          fun = 'mean')
pred.agg = pred.agg[, colnames(w.te.agg)]
pd = acrossRowCor_plotdata(pred = pred.agg, goldstandard = w.te.agg)
# saveRDS(pd, paste0(rdir, 'perf/pred_varmedian_acrosssample_pcc.rds'))
# saveRDS(pd, paste0(rdir, 'perf/pred_allcpg_acrosssample_pcc.rds'))
saveRDS(pd, paste0(rdir, 'perf/', method, '/', key, '_acrosssample_pcc.rds'))

id = complete.cases(pred.agg) & complete.cases(w.te.agg)
pd2 = acrossRowCor_plotdata(pred = t(pred.agg[id,]), 
                            goldstandard = t(w.te.agg[id,]))
# saveRDS(pd2, paste0(rdir, 'perf/pred_varmedian_acrosscpg_pcc.rds'))
saveRDS(pd2, paste0(rdir, 'perf/', method, '/', key, '_acrosscpg_pcc.rds'))

