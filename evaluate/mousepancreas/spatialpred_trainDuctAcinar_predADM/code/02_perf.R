setwd('/home/whou10/data/whou/metpred/')
source('/home/whou10/scratch16/whou10/resource/startup.R')
print(s <- commandArgs(trailingOnly = T)[1]) ## s = 'A72_F_LM_D7'
gs <- readRDS(paste0('evaluate/mousepancreas/goldstandard/res/goldstandard_weightedSum_0.95/', s, '.rds'))
print(str(gs))
dir.create(rdir <- paste0('evaluate/mousepancreas/spatialpred_trainDuctAcinar_predADM/perf/', s, '/'), recursive = T, showWarnings = F)
allf = list.files(paste0('evaluate/mousepancreas/spatialpred_trainDuctAcinar_predADM/res/', s, '/'))
cpgrange = rev(sub('.rds','',sub('pred_cpg_sd_', '', allf)))
cpglist <- readRDS(paste0('evaluate/mousepancreas/celltypespecific_cv/cpg_group/', s, '.rds'))
for (cr in cpgrange){
  print(cr)
  if (!file.exists(paste0(rdir, 'cpg_sd_', cr, '_acrossspot_pcc.rds')) | !file.exists(paste0(rdir, 'cpg_sd_', cr, '_acrosscpg_pcc.rds'))){
    print('Calculating PCC...')
    gs9 = gs[cpglist[[cr]],]
    pred = readRDS(paste0('evaluate/mousepancreas/spatialpred_trainDuctAcinar_predADM/res/', s, '/pred_cpg_sd_', cr, '.rds'))
    pred2 = pred[, colnames(gs)]
    cc <- acrossRowCor_plotdata(pred2, gs9)
    cc.acrossCpG <- acrossRowCor_plotdata(t(pred2), t(gs9))
    saveRDS(cc, paste0(rdir, 'cpg_sd_', cr, '_acrossspot_pcc.rds'))
    saveRDS(cc.acrossCpG, paste0(rdir, 'cpg_sd_', cr, '_acrosscpg_pcc.rds'))  
  }
}


