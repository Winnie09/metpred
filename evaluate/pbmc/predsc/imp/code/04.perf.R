source('/home/whou10/scratch16/whou10/resource/myfunc/01_function.R')
ddir <- '/home/whou10/data/whou/metpred/data/pbmc/imp/res/'
w.te = readRDS(paste0(ddir, 'wgbs_test.rds'))

rdir = '/home/whou10/data/whou/metpred/evaluate/pbmc/predsc/imp/'
pred = readRDS(paste0(rdir, 'res/pred/pred_var_100000_CpG.rds'))
w.te2 = w.te[rownames(pred),]

pred.ct = sub(':.*', '', colnames(pred))
pred.agg = aggregatefunc2(d = pred, by = pred.ct, fun = 'mean')
# num [1:1644, 1:7] 0.816 0.769 0.89 0.298 0.92 ...
# - attr(*, "dimnames")=List of 2
# ..$ : chr [1:1644] "chr22_46144556" "chr10_30482068" "chr18_79422099" "chr2_223918843" ...
# ..$ : chr [1:7] "b_cells" "cd14_monocytes" "cd56_nk" "cytotoxic_t" ...
w.te.agg = aggregatefunc2(d = w.te2, 
                          by = c('cd56_nk', rep('cd14_monocytes', 3), 'b_cells', paste0('t_cells', 1:4)), 
                          fun = 'mean')

pd = acrossRowCor_plotdata(pred = pred.agg, goldstandard = w.te.agg)

saveRDS(pd, paste0(rdir, 'perf/acrosssample_pcc.rds'))

library(ggplot2)
png(paste0(rdir, 'plot/acrosssample_pcc.png'),res = 300,
    width = 650, height = 650)
acrossRowCor_plot(pd[complete.cases(pd),],
                  ylab = 'across-sample PCC')
dev.off()

id = complete.cases(pred.agg) & complete.cases(w.te.agg)
pd2 = acrossRowCor_plotdata(pred = t(pred.agg[id,]), 
                            goldstandard = t(w.te.agg[id,]))
saveRDS(pd2, paste0(rdir, 'perf/acrosscpg_pcc.rds'))

png(paste0(rdir, 'plot/acrosscpg_pcc_diffcpg.png'),res = 300,
    width = 650, height = 650)
acrossRowCor_plot(pd2[complete.cases(pd2),],
                  ylab = 'across-CpG PCC')
dev.off()


