source('/home/whou10/scratch16/whou10/resource/startup.R')
ddir <- '/home/whou10/data/whou/metpred/evaluate/mousepancreas/spatialpred/imp/perf/A46_F_P_K4_D4/'
af = list.files(ddir, 'acrosscpg_pcc_spotset')
res <- lapply(af, function(f){
  tmp = readRDS(paste0(ddir, f))  
})
str(res)  
res <- do.call(rbind, res)
str(res)
pdf('/home/whou10/data/whou/metpred/evaluate/mousepancreas/spatialpred/imp/plot/perf/acrosscpg_pcc_vs_sd_in_gs.pdf', width = 3, height = 2.7)
acrossRowCor_plot(res, ylab = 'across-CpG PCC')
dev.off()

