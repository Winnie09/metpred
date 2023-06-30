gs <- readRDS('evaluate/mousepancreas/goldstandard/res/goldstandard_weightedSum_0.95/A46_F_P_K4_D4.rds')
allf = list.files('evaluate/mousepancreas/spatialpred/imp/res/A46_F_P_K4_D4/')
dir.create(rdir <- 'evaluate/mousepancreas/spatialpred/imp/perf/A46_F_P_K4_D4/', recursive = T)
cpgrange = rev(sub('.rds','',sub('pred_cpg_sd_', '', allf)))
for (cr in cpgrange){
  gs9 = gs[cpglist[[cr]],]
  pred = readRDS(paste0('evaluate/mousepancreas/spatialpred/imp/res/A46_F_P_K4_D4/pred_cpg_sd_', cr, '.rds'))
  pred2 = pred[, colnames(gs)]
  cc <- acrossRowCor_plotdata(pred2, gs9)
  cc.acrossCpG <- acrossRowCor_plotdata(t(pred2), t(gs9))
  saveRDS(cc, paste0(rdir, 'cpg_sd_', cr, '_acrossspot_pcc.rds'))
  saveRDS(cc.acrossCpG, paste0(rdir, 'cpg_sd_', cr, '_acrosscpg_pcc.rds'))
}

