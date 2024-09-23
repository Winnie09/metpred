# key = 'pred_var_10000_CpG'
# key = 'pred_varmedian_0673446_CpG'
key = 'allcpg'
source('/home/whou10/scratch16/whou10/resource/myfunc/01_function.R')
ddir <- '/home/whou10/data/whou/metpred/data/pbmc/imp/res/'
rdir = '/home/whou10/data/whou/metpred/evaluate/pbmc/predsc/imp/'
pred.agg <- readRDS(paste0(rdir, 'eval_pred_ctmean.rds'))
w.te.agg <- readRDS(paste0(rdir, 'eval_goldstandard_ctmean.rds'))

pd = acrossRowCor_plotdata(pred = pred.agg, goldstandard = w.te.agg)
saveRDS(pd, paste0(rdir, 'perf/', key, '_acrosssample_pcc.rds'))

library(ggplot2)
# png(paste0(rdir, 'plot/pred_varmedian_acrosssample_pcc.png'),res = 300,
#     width = 650, height = 650)
png(paste0(rdir, 'plot/', key, '_acrosssample_pcc.png'),res = 300,
    width = 650, height = 650)
acrossRowCor_plot(pd[complete.cases(pd),],
                  ylab = 'across-sample PCC',
                  xlab = 'testing measured value sd')
dev.off()

## merge intervals that have less then 5e4 CpGs, i.e. merge (0.35 to 0.1] 
levels(pd[,1])[levels(pd[,1]) %in% paste0('(',seq(0.35,0.95,0.05), ',',seq(0.35,0.95,0.05)+0.05,']')] <- '(0.35,1]'

library(ggplot2)
png(paste0(rdir, 'plot/', key, '_acrosssample_pcc_merge_interval.png'),res = 300,
    width = 650, height = 650)
acrossRowCor_plot(pd[complete.cases(pd),],
                  ylab = 'across-sample PCC',
                  xlab = 'testing measured value sd')
dev.off()



##
pd2 = acrossRowCor_plotdata(pred = t(pred.agg), goldstandard = t(w.te.agg))
str(pd2)
saveRDS(pd2, paste0(rdir, 'perf/', key, '_acrossCpG_pcc.rds'))

library(ggplot2)
png(paste0(rdir, 'plot/', key, '_acrossCpG_pcc.png'),res = 300,
    width = 650, height = 650)
acrossRowCor_plot(pd2[complete.cases(pd2),],
                  ylab = 'across-CpG PCC',
                  xlab = 'testing measured value sd',
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
                  ylab = 'across-CpG PCC',
                  xlab = 'testing measured value sd',
                  point = T,
                  point.size = 2,
                  point.stroke = 0,
                  point.alpha = 0.5,
                  point.width = 0.1)
dev.off()


pd[100:110,]
par(mfrow=c(1,2))
cpg = 'chr16_22435056'
plot(pred.agg[cpg,], w.te.agg[cpg,], xlab='pred', ylab='measured', pch = 20, main=cpg)

cpg = 'chr16_22435145'
plot(pred.agg[cpg,], w.te.agg[cpg,], xlab='pred', ylab='measured', pch = 20, main=cpg)

identical(rownames(pd), rownames(pred.agg))
smoothScatter(pd[,2], pred.agg[,1], ylab='pred DNAm', xlab='across-sample PCC', pch = 20, size = 0.001)

id = complete.cases(pred.agg) & complete.cases(w.te.agg)
pd2 = acrossRowCor_plotdata(pred = t(pred.agg[id,]), 
                            goldstandard = t(w.te.agg[id,]))
# saveRDS(pd2, paste0(rdir, 'perf/pred_varmedian_acrosscpg_pcc.rds'))
saveRDS(pd2, paste0(rdir, 'perf/', key, '_acrosscpg_pcc.rds'))

# png(paste0(rdir, 'plot/pred_varmedian_acrosscpg_pcc_diffcpg.png'),res = 300,
#     width = 650, height = 650)

png(paste0(rdir, 'plot/', key, '_acrosscpg_pcc_diffcpg.png'),res = 300,
    width = 650, height = 650)
acrossRowCor_plot(pd2[complete.cases(pd2),],
                  ylab = 'across-CpG PCC',
                  xlab = 'testing measured value sd')
dev.off()


##########################



