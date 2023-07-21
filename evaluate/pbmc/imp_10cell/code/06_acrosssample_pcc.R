source('/home/whou10/scratch16/whou10/resource/myfunc/01_function.R')
ddir <- '/home/whou10/data/whou/metpred/data/pbmc/imp_10cell/res/'
w.te = readRDS(paste0(ddir, 'wgbs_test.rds'))
pd_acrosssample <- list()
rdir1 = '/home/whou10/data/whou/metpred/evaluate/pbmc/imp_10cell/res/predall/'
for (i in  list.files(rdir1)){
  print(i)
  pred = readRDS(paste0(rdir1, i))
  str(pred)
  w.te2 = w.te[rownames(pred),]
  str(w.te2)
  colnames(w.te2)
  pred.ct = sub(':.*', '', colnames(pred))
  pred.agg = aggregatefunc2(d = pred, by = pred.ct, fun = 'mean')
  w.te.agg = aggregatefunc2(d = w.te2, 
                            by = c('cd56_nk', rep('cd14_monocytes', 3), 'b_cells', paste0('t_cells', 1:4)), 
                            fun = 'mean')
  pd_acrosssample[[sub('.rds','',i)]] <- acrossRowCor_plotdata(pred = pred.agg, goldstandard = w.te.agg)
}
pd = do.call(rbind, pd_acrosssample)  
saveRDS(pd, '/home/whou10/data/whou/metpred/evaluate/pbmc/imp_10cell/perf/acrosssample_pcc.rds')
saveRDS(pd_acrosssample, '/home/whou10/data/whou/metpred/evaluate/pbmc/imp_10cell/perf/acrosssample_pcc_list.rds')

png('/home/whou10/data/whou/metpred/evaluate/pbmc/imp_10cell/plot/acrosssample_pcc.png',res = 300,
    width = 650, height = 650)
acrossRowCor_plot(pd[complete.cases(pd),],
                  ylab = 'across-sample PCC')
dev.off()

pd_acrosssample2 = list()
for (i in 1:length(pd_acrosssample)){
  m = pd_acrosssample[[i]]
  m[,1] = paste0('(', names(pd_acrosssample)[i], ']')
  pd_acrosssample2[[i]] = m
}
pd2 = do.call(rbind, pd_acrosssample2)  
png('/home/whou10/data/whou/metpred/evaluate/pbmc/imp_10cell/plot/acrosssample_pcc_training_sd.png',res = 300,
    width = 650, height = 650)
acrossRowCor_plot(pd2[complete.cases(pd2),],
                  ylab = 'across-sample PCC',
                  xlab = 'Training measured value sd')
dev.off()

