source('/home/whou10/scratch16/whou10/resource/myfunc/01_function.R')
ddir <- '/home/whou10/data/whou/metpred/data/pbmc/res/'
w.te = readRDS(paste0(ddir, 'wgbs_test.rds'))
w.te2 = w.te[rownames(pred),]
pd = acrossRowCor_plotdata(pred = pred, goldstandard = w.te2)
saveRDS(pd, '/home/whou10/data/whou/metpred/evaluate/pbmc/perf/acrosssample_pcc.rds')

png('/home/whou10/data/whou/metpred/evaluate/pbmc/plot/acrosssample_pcc.png',res = 300,
    width = 650, height = 650)
acrossRowCor_plot(pd[complete.cases(pd),],
                  ylab = 'across-sample PCC')
dev.off()

id = complete.cases(pred) & complete.cases(w.te2)
pd2 = acrossRowCor_plotdata(pred = t(pred[id,]), 
                            goldstandard = t(w.te2[id,]))
saveRDS(pd2, '/home/whou10/data/whou/metpred/evaluate/pbmc/perf/acrosscpg_pcc.rds')

png('/home/whou10/data/whou/metpred/evaluate/pbmc/plot/acrosscpg_pcc.png',res = 300,
    width = 650, height = 650)
acrossRowCor_plot(pd2[complete.cases(pd2),],
                  ylab = 'across-CpG PCC')
dev.off()

