setwd('/home/whou10/data/whou/metpred/')
ddir = 'evaluate/mousepancreas/spatialpred/imp/perf/A46_F_P_K4_D4/'
pdir = 'evaluate/mousepancreas/spatialpred/imp/plot/perf/'
allf = list.files(ddir)
allf = allf[grep('acrossspot_pcc', allf)]
cclist <- lapply(allf, function(f){
  tmp <- readRDS(paste0(ddir, f))
})
cc <- do.call(rbind, cclist)  
pdf(paste0(pdir, 'acrossspot_pcc_vs_cpg_sd_in_gs.pdf'),width = 3, height = 2.7)
acrossRowCor_plot(plotdata = cc, 
                  ylab = 'across-spot PCC',
                  xlab = 'CpG sd in goldstandard')
dev.off()
###
allf = list.files(ddir)
allf = allf[grep('acrossspot_pcc', allf)]
cclist <- lapply(allf, function(f){
  tmp <- readRDS(paste0(ddir, f))
  tmp[,1] <- paste0('(','',sub('_acrossspot_pcc.rds','',sub('cpg_sd_', '', f)), ']')
  tmp
})
cc <- do.call(rbind, cclist)  
pdf(paste0(pdir, 'acrossspot_pcc_vs_cpg_sd_in_training.pdf'), width = 3, height = 2.7)
acrossRowCor_plot(plotdata = cc, 
                  ylab = 'across-spot PCC',
                  xlab = 'CpG sd in training samples')
dev.off()

