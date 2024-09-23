# key = 'pred_var_100000_CpG'
# key = 'pred_varmedian_0673446_CpG'
key = 'allcpg'
source('/home/whou10/scratch16/whou10/resource/myfunc/01_function.R')
ddir <- '/home/whou10/data/whou/metpred/evaluate/pbmc/predsc/imp/res/eval/'
rdir = '/home/whou10/data/whou/metpred/evaluate/pbmc/predsc/imp/'
pred.agg <- readRDS(paste0(ddir, 'eval_pred_ctmean.rds'))
w.te.agg <- readRDS(paste0(ddir, 'eval_goldstandard_ctmean.rds'))

pd = acrossRowCor_plotdata(pred = pred.agg, goldstandard = w.te.agg)
saveRDS(pd, paste0(rdir, 'perf/ramp/', key, '_acrosssample_pcc.rds'))

library(ggplot2)
png(paste0(rdir, 'plot/', key, '_acrosssample_pcc.png'),res = 300,
    width = 650, height = 650)
acrossRowCor_plot(pd[complete.cases(pd),],
                  ylab = 'Across-sample PCC',
                  xlab = 'Across-sample variability')
dev.off()

## merge intervals that have less then 5e4 CpGs, i.e. merge (0.35 to 0.1] 
levels(pd[,1])[levels(pd[,1]) %in% paste0('(',seq(0.35,0.95,0.05), ',',seq(0.35,0.95,0.05)+0.05,']')] <- '(0.35,1]'

library(ggplot2)
png(paste0(rdir, 'plot/', key, '_acrosssample_pcc_merge_interval.png'),res = 300,
    width = 650, height = 650)
acrossRowCor_plot(pd[complete.cases(pd),],
                  ylab = 'Across-sample PCC',
                  xlab = 'Across-sample variability')
dev.off()


## acrosscpg pcc
pd2 = acrossRowCor_plotdata(pred = t(pred.agg), goldstandard = t(w.te.agg))
saveRDS(pd2, paste0(rdir, 'perf/ramp/', key, '_acrosscpg_pcc.rds'))

library(ggplot2)
png(paste0(rdir, 'plot/', key, '_acrossCpG_pcc.png'),res = 300,
    width = 650, height = 650)
acrossRowCor_plot(pd2[complete.cases(pd2),],
                  ylab = 'Across-CpG PCC',
                  xlab = 'Across-CpG variability',
                  point = T,
                  point.size = 2,
                  point.stroke = 0,
                  point.alpha = 0.5,
                  point.width = 0.1)
dev.off()


## merge all intervals
pd2[,1] = as.character(pd2[,1])
pd2[,1] <- '(0.2,0.3]'
png(paste0(rdir, 'plot/', key, '_acrossCpG_pcc_merge_interval.png'),res = 300,
    width = 650, height = 650)
acrossRowCor_plot(pd2[complete.cases(pd2),],
                  ylab = 'Across-CpG PCC',
                  xlab = 'Across-CpG variability',
                  point = T,
                  point.size = 2,
                  point.stroke = 0,
                  point.alpha = 0.5,
                  point.width = 0.1)
dev.off()
