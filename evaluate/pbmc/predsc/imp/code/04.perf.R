key = 'pred_var_10000_CpG'
# key = 'pred_varmedian_0673446_CpG'

acrossRowCor_plot <- function(plotdata = NA, 
                              xlab = 'Testing measured value sd', 
                              ylab = 'across-row PCC', 
                              title = '', 
                              savefile = FALSE,
                              filename = 'acrossRow_pcc.pdf',
                              width = 2.5,
                              facet.var = NA,
                              facet.nrow = 3){
  ## plotdata: dataframe. first column is standard deviation (x-axis), second column is
  ## correlation (y-axis). 
  p <- ggplot(data = plotdata, aes(y=cor,x=sd,fill=sd)) + 
    geom_violin(alpha=0.2, scale = 'width') + 
    geom_boxplot(alpha=0.3, width = 0.2, outlier.size = 0.01) + 
    # geom_jitter(alpha=0.1, width = 0.1, stroke = 0, size = 0.1, color = 'grey3') + 
    theme_classic() + 
    scale_fill_manual(values=rainbow(length(unique(plotdata$sd)))) + 
    xlab(xlab) + 
    ylab(ylab) + 
    theme(legend.position = 'none', axis.text.x = element_text(angle = 45, hjust = 1)) + 
    ggtitle(title)
  if (!is.na(facet.var)){
    eval(parse(text=paste0('p <- p + facet_wrap(~',facet.var,',nrow=', facet.nrow,')')))
  }
  if (savefile) {
    pdf(filename,width=width,height=2.5)
    print(p)
    dev.off()
  } else {
    return(p)
  }  
}

source('/home/whou10/scratch16/whou10/resource/myfunc/01_function.R')
ddir <- '/home/whou10/data/whou/metpred/data/pbmc/imp/res/'
w.te = readRDS(paste0(ddir, 'wgbs_test.rds'))

rdir = '/home/whou10/data/whou/metpred/evaluate/pbmc/predsc/imp/'
# pred = readRDS(paste0(rdir, 'res/pred/pred_varmedian_30673446_CpG.rds'))
# pred = readRDS(paste0(rdir, 'res/pred/pred_allcpg.rds'))
pred = readRDS(paste0(rdir, 'res/pred/pred_var_100000_CpG.rds'))
w.te2 = w.te[rownames(pred),]
pred.ct = sub(':.*', '', colnames(pred))
pred.ct[grepl('_t', pred.ct)] = 't_cells'
pred.agg = aggregatefunc2(d = pred, by = pred.ct, fun = 'mean')
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
saveRDS(pd, paste0(rdir, 'perf/', key, '_acrosssample_pcc.rds'))

library(ggplot2)
# png(paste0(rdir, 'plot/pred_varmedian_acrosssample_pcc.png'),res = 300,
#     width = 650, height = 650)
png(paste0(rdir, 'plot/', key, '_acrosssample_pcc.png'),res = 300,
    width = 650, height = 650)
acrossRowCor_plot(pd[complete.cases(pd),],
                  ylab = 'across-sample PCC')
dev.off()




head(pd)
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
                  ylab = 'across-CpG PCC')
dev.off()



