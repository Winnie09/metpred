source('/home/whou10/scratch16/whou10/resource/startup.R')

ddir0 <- 'evaluate/mousepancreas/spatialpred/imp/perf/'
alls = list.files(ddir0)
s = alls[1]
res2 = lapply(alls, function(s){
  print(s)
  ddir = paste0(ddir0, s, '/')
  af = list.files(ddir, 'acrosscpg_pcc_spotset')
  res <- lapply(af, function(f){
    tmp = readRDS(paste0(ddir, f))  
    tmp2 = data.frame(tmp, sample = s)
  })
  res <- do.call(rbind, res)
})
str(res2) 
res3 = do.call(rbind, res2)
str(res3)

p <- acrossRowCor_plot(res3, ylab = 'across-CpG PCC', facet.var = 'sample')
ggsave('/home/whou10/data/whou/metpred/evaluate/mousepancreas/spatialpred/imp/plot/perf/all/acrosscpg_pcc_vs_sd_in_gs_all.png', 
       p,
       width = 5.6,
       height = 4.5,
       dpi = 300)

sels <- c('A46_F_P_K4_D4', 'A47_F_P_K4_D7', 'A72_F_LM_D7', 'A73_F_LM_D4')
p <- acrossRowCor_plot(res3[res3[,3] %in% sels, ], 
                       ylab = 'across-CpG PCC', 
                       facet.var = 'sample',
                       facet.nrow = 1)
ggsave('/home/whou10/data/whou/metpred/evaluate/mousepancreas/spatialpred/imp/plot/perf/all/acrosscpg_pcc_vs_sd_in_gs_select.png', 
       p,
       width = 5.6,
       height = 2,
       dpi = 300)


