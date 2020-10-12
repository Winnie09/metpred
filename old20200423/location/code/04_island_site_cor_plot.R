setwd('/home-4/whou10@jhu.edu/scratch/Wenpin/')
# setwd('/Users/wenpinhou/Dropbox/')
library(ggplot2)
key = as.character(commandArgs(trailingOnly = T)[1]) ## wgbs, array
# key = 'wgbs'
chrn = as.character(commandArgs(trailingOnly = T)[2])
# chrn = 'chr21'
df = readRDS(paste0('./metpred/location/plot/',key,'/island_site_cor_',chrn,'_pd.rds'))
df$type = as.character(df$type)
df$type = ifelse(df$type=='two sites(same island)','two sites(same CGI)', 'two sites(different CGI)')
df$type = factor(df$type, levels=c('two sites(same CGI)', 'two sites(different CGI)'))
head(df)
v = df[df$type == unique(df$type)[1],1]
u = df[df$type == unique(df$type)[2],1]
pval = t.test(v,u, alternative = "two.sided")$p.value
pdf(paste0('./metpred/location/plot/',key,'/island_site_cor_',chrn,'.pdf'),width=4,heigh=2)
ggplot(data=df,aes(x=type, y=correlation,fill=type)) + geom_violin(alpha=0.2,trim=T)  + 
  geom_boxplot(width=0.1,outlier.shape = NA) +
  ggtitle(ifelse((pval-1e7)<0,expression('p.value<'~10^7),paste0('p.value=',formatC(pval, format = "e", digits = 3))))+
  theme_classic() + xlab(paste0('Chromosome ',sub('chr','',chrn)))+theme(legend.position = 'none', axis.text.x=element_text(color='black',size=11), axis.title = element_text(size=11)) + ylab('Spearman correlation')
dev.off()


