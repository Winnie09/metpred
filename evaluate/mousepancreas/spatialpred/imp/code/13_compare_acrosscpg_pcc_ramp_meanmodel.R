########################################
## merge results
ddir = '/home/whou10/data/whou/metpred/evaluate/mousepancreas/spatialpred/imp/perf/meanmodel/'
alls = list.files(ddir)
s = alls[2]
for (s in alls){
  print(s)
  allf = list.files(paste0(ddir, s), pattern = 'spotset')
  if (length(allf) > 0) {
    d <- sapply(allf, function(f){
      tmp = readRDS(paste0(ddir, s, '/', f))
    }, simplify = F)
    d = unlist(d)
    names(d) = sub('.*.rds.', '', names(d))
  }
  str(d)
  saveRDS(d, paste0(ddir, s, '/acrosscpg_pcc_allspot.rds'))
}


#########################################
## plot results
source('/home/whou10/scratch16/whou10/resource/startup.R')

ddir0 <- 'evaluate/mousepancreas/spatialpred/imp/perf/ramp/'
ddir.mean = '/home/whou10/data/whou/metpred/evaluate/mousepancreas/spatialpred/imp/perf/meanmodel/'
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
  res.m = readRDS(paste0(ddir.mean, s, '/acrosscpg_pcc_allspot.rds'))    
  res.m = res.m[rownames(res)]
  res.tmp = res
  res.tmp[,2] = res.m
  res.all = rbind(data.frame(res, model = 'ramp'), data.frame(res.tmp, model = 'meanmodel'))
})
str(res2) 
res3 = do.call(rbind, res2)
str(res3)

## use bar plot
library(scales)
result <- aggregate(cor ~ sample + model, data = res3[res3[,3] %in% sels, ], FUN = median)
p <- ggplot()+
  geom_bar(data = result, aes(y=cor, x=model, fill=model, color=model), stat='identity', alpha = 0.1) +
  geom_jitter(data = res3[res3[,3] %in% sels,], aes(y=cor,x=model,fill=model, color=model),
              size = .3, stroke = 0, alpha = 0.5, width = 0.3) + 
  theme_classic() + 
  xlab('Method') + 
  ylab('Across-CpG PCC') + 
  theme(legend.position = 'none', axis.text.x = element_text(angle = 45, hjust = 1)) + 
  ggtitle('') +
  facet_wrap(~sample, nrow = 1) +
  scale_y_continuous(limits=c(0.9,1),oob = rescale_none) +
  scale_fill_brewer(palette = 'Set1') +
  scale_color_brewer(palette = 'Set1')

ggsave('/home/whou10/data/whou/metpred/evaluate/mousepancreas/spatialpred/imp/plot/perf/all/acrosscpg_pcc_vs_sd_in_gs_select.png', 
       p,
       width = 5.6,
       height = 2.4,
       dpi = 300)

###################
## use boxplot
p <- ggplot(data = res3, aes(y=cor,x=sd,fill=model)) + 
  geom_violin(alpha=0.2, scale = 'width') + 
  geom_boxplot(alpha=0.3, width = 0.2, outlier.size = 0.01) + 
  theme_classic() + 
  scale_fill_manual(values=rainbow(length(unique(res3$sd)))) + 
  xlab('testing measured value sd') + 
  ylab('across-CpG PCC') + 
  theme(legend.position = 'none', axis.text.x = element_text(angle = 45, hjust = 1)) + 
  ggtitle('') +
  facet_wrap(~sample)
  
ggsave('/home/whou10/data/whou/metpred/evaluate/mousepancreas/spatialpred/imp/plot/perf/all/acrosscpg_pcc_vs_sd_in_gs_all_boxplot.png', 
       p,
       width = 5.6,
       height = 4.5,
       dpi = 300)

sels <- c('A46_F_P_K4_D4', 'A47_F_P_K4_D7', 'A72_F_LM_D7', 'A73_F_LM_D4')


p <- ggplot(data = res3[res3[,3] %in% sels, ], aes(y=cor,x=sd,fill=model, color=model)) + 
  #geom_violin(alpha=0.2, scale = 'width') + 
  geom_boxplot(alpha=0.3, width = 0.2, outlier.size = 0.01) + 
  theme_classic() + 
  #scale_fill_manual(values=rainbow(length(unique(res3$sd)))) + 
  xlab('testing measured value sd') + 
  ylab('across-CpG PCC') + 
  theme(legend.position = 'none', axis.text.x = element_text(angle = 45, hjust = 1)) + 
  ggtitle('') +
  facet_wrap(~sample, nrow = 1)

ggsave('/home/whou10/data/whou/metpred/evaluate/mousepancreas/spatialpred/imp/plot/perf/all/acrosscpg_pcc_vs_sd_in_gs_all_boxplot.png', 
       p,
       width = 5.6,
       height = 2.4,
       dpi = 300)

