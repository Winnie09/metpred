##################
rm(list=ls())
source('/home/whou10/scratch16/whou10/resource/startup.R')
rdir <- '/home/whou10/data/whou/metpred/evaluate/mousepancreas/celltypespecific_cv/perf/all/'
df = data.frame(method = c('ramp', 'neargene', 'permu'), fullname = c('Ramp', 'Nearest gene', 'Permute'), stringsAsFactors = FALSE)
m = lapply(1:nrow(df), function(i){
  tmp = data.frame(readRDS(paste0(rdir, df[i,1], '_across_sample_mse.rds')), method = df[i,2])
})
mse = do.call(rbind, m)

colv = c('red', 'orange', 'royalblue')
pal <- c('Ramp'='red','Nearest gene'='orange','Permute'='royalblue')
pdf('/home/whou10/data/whou/metpred/evaluate/mousepancreas/celltypespecific_cv/plot/allmethods_across_sample_mse.pdf',
    width = 4, height = 1.5)
ggplot(data = mse, aes(y=mse,x=sd,fill=method)) + 
  geom_boxplot(alpha=0.8, width = 0.6, outlier.size = 0) + 
  scale_fill_manual(values=pal) + 
  # scale_fill_manual(values=rainbow(length(unique(plotdata$sd)))) +
  xlab('Across-sample CpG variability') +
  ylab('Across-sample MSE') +
  theme(axis.text.x = element_text(angle = 45, hjust = 0.5)) 
dev.off()


