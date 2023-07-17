source('/home/whou10/scratch16/whou10/resource/myfunc/01_function.R')
ddir <- '/home/whou10/data/whou/metpred/data/pbmc/imp/res/'
w.te = readRDS(paste0(ddir, 'wgbs_test.rds'))
pred = readRDS('/home/whou10/data/whou/metpred/evaluate/pbmc/res/pred/pred_top1e4_mono_t_diff_greater0.2_tr_te.rds')
w.te2 = w.te[rownames(pred),]

pred.ct = sub(':.*', '', colnames(pred))
pred.agg = aggregatefunc2(d = pred, by = pred.ct, fun = 'mean')
w.te.agg = aggregatefunc2(d = w.te2, 
                          by = c('cd56_nk', rep('cd14_monocytes', 3), 'b_cells', paste0('t_cells', 1:4)), 
                          fun = 'mean')

pd = acrossRowCor_plotdata(pred = pred.agg, goldstandard = w.te.agg)
saveRDS(pd, '/home/whou10/data/whou/metpred/evaluate/pbmc/perf/acrosssample_pcc.rds')

png('/home/whou10/data/whou/metpred/evaluate/pbmc/plot/acrosssample_pcc.png',res = 300,
    width = 650, height = 650)
acrossRowCor_plot(pd[complete.cases(pd),],
                  ylab = 'across-sample PCC')
dev.off()

id = complete.cases(pred.agg) & complete.cases(w.te.agg)
pd2 = acrossRowCor_plotdata(pred = t(pred.agg[id,]), 
                            goldstandard = t(w.te.agg[id,]))
saveRDS(pd2, '/home/whou10/data/whou/metpred/evaluate/pbmc/perf/acrosscpg_pcc.rds')

png('/home/whou10/data/whou/metpred/evaluate/pbmc/plot/acrosscpg_pcc.png',res = 300,
    width = 650, height = 650)
acrossRowCor_plot(pd2[complete.cases(pd2),],
                  ylab = 'across-CpG PCC')
dev.off()


