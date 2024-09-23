## performance
# key = 'pred_var_10000_CpG'
# key = 'pred_varmedian_0673446_CpG'
key = 'allcpg'
source('/home/whou10/scratch16/whou10/resource/myfunc/01_function.R')
ddir <- '/home/whou10/data/whou/metpred/data/pbmc/imp/res/'
rdir = '/home/whou10/data/whou/metpred/evaluate/pbmc/predsc/imp/res/eval/'
w.te = readRDS(paste0(ddir, 'wgbs_test.rds'))
pred.agg = readRDS(paste0('/home/whou10/data/whou/metpred/evaluate/pbmc/predsc/imp/res/pb/', 'pred_allcpg.rds'))
w.te2 = w.te[rownames(pred.agg),]
w.te.agg = aggregatefunc2(d = w.te2, 
                          by = c('cd56_nk', rep('cd14_monocytes', 3), 'b_cells', rep('t_cells', 4)), 
                          fun = 'mean')
pred.agg = pred.agg[, colnames(w.te.agg)]
saveRDS(pred.agg, paste0(rdir, 'eval_pred_ctmean.rds'))
saveRDS(w.te.agg, paste0(rdir, 'eval_goldstandard_ctmean.rds'))


## ==========
key = 'pred_var_10000_CpG'
source('/home/whou10/scratch16/whou10/resource/myfunc/01_function.R')
ddir <- '/home/whou10/data/whou/metpred/data/pbmc/imp/res/'
rdir = '/home/whou10/data/whou/metpred/evaluate/pbmc/predsc/imp/res/eval/'
w.te = readRDS(paste0(ddir, 'wgbs_test.rds'))
pred <- readRDS('/home/whou10/data/whou/metpred/evaluate/pbmc/predsc/imp/res/pred/pred_var_100000_CpG.rds')
w.te2 = w.te[rownames(pred),]
w.te3 = aggregatefunc2(d = w.te2, 
                          by = c('cd56_nk', rep('cd14_monocytes', 3), 'b_cells', rep('t_cells', 4)), 
                          fun = 'mean')
pred = pred[, colnames(w.te.agg)]
saveRDS(pred, paste0(rdir, 'eval_pred_ctmean_var_100000CpG.rds'))
saveRDS(w.te3, paste0(rdir, 'eval_goldstandard_ctmean_var_100000CpG.rds'))

