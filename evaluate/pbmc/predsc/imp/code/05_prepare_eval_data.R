rm(list=ls())
## prepare evaluation data: goldstandard, and predicted sc pooled pb
## for all cpg prediction
source('/home/whou10/scratch16/whou10/resource/myfunc/01_function.R')
ddir <- '/home/whou10/data/whou/metpred/data/pbmc/imp/res/'
rdir = '/home/whou10/data/whou/metpred/evaluate/pbmc/predsc/imp/res/eval/'
w.te = readRDS(paste0(ddir, 'wgbs_test.rds'))
pred.agg = readRDS('/home/whou10/data/whou/metpred/evaluate/pbmc/predsc/imp/res/pb/pred_allcpg.rds')
w.te2 = w.te[rownames(pred.agg),]
w.te.agg = aggregatefunc2(d = w.te2, 
                          by = c('cd14_monocytes', 'b_cells', rep('t_cells', 3)), 
                          fun = 'mean')
pred.agg = pred.agg[, colnames(w.te.agg)]
saveRDS(pred.agg, paste0(rdir, 'eval_pred_ctmean.rds'))
saveRDS(w.te.agg, paste0(rdir, 'eval_goldstandard_ctmean.rds'))

pdir = '/home/whou10/data/whou/metpred/evaluate/pbmc/predsc/imp/plot/'
png(paste0(pdir, 'compare_pred_and_gs_smoothscatter.png'), width = 2e3, height = 2e3, res = 200)
par(mfrow=c(3,2))
mdf = data.frame(v1 = c(1,2,3,4,5), v2=c(1,2,3,3,3))
for (i in 1:nrow(mdf)){
  print(i)
  smoothScatter(w.te2[,mdf[i,1]], 
                pred.agg[,mdf[i,2]],
                xlab=colnames(w.te2)[mdf[i,1]], 
                ylab=colnames(pred.agg)[mdf[i,2]], 
                main = round(cor(w.te2[,mdf[i,1]], pred.agg[,mdf[i,2]]), 3))
  
}
dev.off()



## ====================================================================
## for top 10e5 cpg with highest var
source('/home/whou10/scratch16/whou10/resource/myfunc/01_function.R')
ddir <- '/home/whou10/data/whou/metpred/data/pbmc/imp/res/'
rdir = '/home/whou10/data/whou/metpred/evaluate/pbmc/predsc/imp/res/eval/'
w.te = readRDS(paste0(ddir, 'wgbs_test.rds'))
pred <- readRDS('/home/whou10/data/whou/metpred/evaluate/pbmc/predsc/imp/res/pred/pred_var_100000_CpG.rds')
pred.ct = sub(':.*', '', colnames(pred))
pred.ct[grepl('_t', pred.ct)] = 't_cells'
pred.agg = aggregatefunc2(d = pred, by = pred.ct, fun = 'mean')

w.te2 = w.te[rownames(pred),]
w.te.agg = aggregatefunc2(d = w.te2, 
                          by = c('cd14_monocytes', 'b_cells', rep('t_cells', 3)), 
                          fun = 'mean')
pred.agg = pred.agg[, colnames(w.te.agg)]
saveRDS(pred.agg, paste0(rdir, 'eval_pred_ctmean_var_100000CpG.rds'))
saveRDS(w.te.agg, paste0(rdir, 'eval_goldstandard_ctmean_var_100000CpG.rds'))

pdir = '/home/whou10/data/whou/metpred/evaluate/pbmc/predsc/imp/plot/'
png(paste0(pdir, 'compare_pred_and_gs_var_100000cpg_smoothscatter.png'), width = 2e3, height = 2e3, res = 200)
par(mfrow=c(3,2))
mdf = data.frame(v1 = c(1,2,3,4,5), v2=c(2,1,4,4,4))
for (i in 1:nrow(mdf)){
  
  smoothScatter(w.te2[,mdf[i,1]], 
                pred.agg[,mdf[i,2]],
                xlab=colnames(w.te2)[mdf[i,1]], 
                ylab=colnames(pred.agg)[mdf[i,2]], 
                main = round(cor(w.te2[,mdf[i,1]], pred.agg[,mdf[i,2]]), 3))
  
}
dev.off()
