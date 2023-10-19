setwd('/home/whou10/data/whou/metpred/')
ddir = 'evaluate/mousepancreas/spatialpred_trainDuctAcinar_predADM/perf/'
alls = list.files(ddir)
source('/home/whou10/scratch16/whou10/resource/startup.R')
for (s in alls){
  print(s)
  ddir = 'evaluate/mousepancreas/spatialpred_trainDuctAcinar_predADM/perf/'
  ddir = paste0(ddir, s, '/')
  dir.create(pdir <- paste0('evaluate/mousepancreas/spatialpred_trainDuctAcinar_predADM/plot/perf/', s, '/'), recursive = T)
  allf = list.files(ddir)
  allf = allf[grep('acrossspot_pcc', allf)]
  cclist <- lapply(allf, function(f){
    tmp <- readRDS(paste0(ddir, f))
  })
  cc <- do.call(rbind, cclist)  
  saveRDS(cc, paste0(ddir, 'acrossspotPCC_vs_cpg_sd_in_gs.rds'))
  pdf(paste0(pdir, 'acrossspot_pcc_vs_cpg_sd_in_gs.pdf'),width = 3, height = 2.7)
  print(acrossRowCor_plot(plotdata = cc, 
                    ylab = 'across-spot PCC',
                    xlab = 'CpG sd in goldstandard',
                    title = s))
  dev.off()
  try(dev.off(dev.list()["RStudioGD"]), silent=TRUE)
  ###
  allf = list.files(ddir)
  allf = allf[grep('acrossspot_pcc', allf)]
  cclist <- lapply(allf, function(f){
    tmp <- readRDS(paste0(ddir, f))
    tmp[,1] <- paste0('(','',sub('_acrossspot_pcc.rds','',sub('cpg_sd_', '', f)), ']')
    tmp
  })
  cc <- do.call(rbind, cclist)  
  saveRDS(cc, paste0(ddir, 'acrossspotPCC_vs_cpg_sd_in_training.rds'))
  pdf(paste0(pdir, 'acrossspot_pcc_vs_cpg_sd_in_training.pdf'), width = 3, height = 2.7)
  print(acrossRowCor_plot(plotdata = cc, 
                    ylab = 'across-spot PCC',
                    xlab = 'CpG sd in training samples',
                    title = s))
  dev.off()
  try(dev.off(dev.list()["RStudioGD"]), silent=TRUE)
}



