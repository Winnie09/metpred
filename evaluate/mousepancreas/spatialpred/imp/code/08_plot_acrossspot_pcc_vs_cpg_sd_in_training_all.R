setwd('/home/whou10/data/whou/metpred/')
ddir = 'evaluate/mousepancreas/spatialpred/imp/perf/ramp/'
alls = list.files(ddir)
# ## only plot four samples
# alls <- sels <- c('A46_F_P_K4_D4', 'A47_F_P_K4_D7', 'A72_F_LM_D7', 'A73_F_LM_D4')

## plot
ccall = lapply(alls, function(s){
  print(s)
  ddir = paste0('evaluate/mousepancreas/spatialpred/imp/perf/ramp/', s, '/')
  cc = readRDS(paste0(ddir, 'acrossspotPCC_vs_cpg_sd_in_training.rds'))
  set.seed(1)
  cc = cc[sample(1:nrow(cc), 1e4), ]
  tmp = data.frame(cc, sample = s)
})  
str(ccall[[1]])  
cc = do.call(rbind, ccall)
source('/home/whou10/scratch16/whou10/resource/startup.R')
p <- ggplot(data = cc, aes(y=cor,x=sd,fill=sd)) + 
  geom_violin(alpha=0.2, scale = 'width') + 
  geom_boxplot(alpha=0.3, width = 0.2, outlier.size = 0.01) + 
  theme_classic() + 
  scale_fill_manual(values=rainbow(length(unique(cc$sd)))) + 
  ylab('across-spot PCC') + 
  xlab('CpG sd in training samples') + 
  theme(legend.position = 'none', axis.text.x = element_text(angle = 45, hjust = 1)) + 
  facet_wrap(~sample)
ggsave('evaluate/mousepancreas/spatialpred/imp/plot/perf/all/acrossspot_pcc_vs_cpg_sd_in_training.png',
       p,
       width = 5.6,
       height = 4.5,
       dpi = 300)

sels <- c('A46_F_P_K4_D4', 'A47_F_P_K4_D7', 'A72_F_LM_D7', 'A73_F_LM_D4')
p <- ggplot(data = cc[cc[,3]%in%sels, ], aes(y=cor,x=sd,fill=sd)) + 
  geom_violin(alpha=0.2, scale = 'width') + 
  geom_boxplot(alpha=0.3, width = 0.2, outlier.size = 0.01) + 
  theme_classic() + 
  scale_fill_manual(values=rainbow(length(unique(cc$sd)))) + 
  ylab('Across-spot PCC') + 
  xlab('CpG sd in training samples') + 
  theme(legend.position = 'none', axis.text.x = element_text(angle = 45, hjust = 1)) + 
  facet_wrap(~sample, nrow = 1) + 
  geom_hline(yintercept=0, color = 'red', linetype = 'dashed')
p
 ggsave('evaluate/mousepancreas/spatialpred/imp/plot/perf/all/acrossspot_pcc_vs_cpg_sd_in_training_sel.png',
        p,
        width = 5.6,
        height = 2,
        dpi = 300)

