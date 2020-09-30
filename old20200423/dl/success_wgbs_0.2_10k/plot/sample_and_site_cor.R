setwd('/Users/wenpinhou/Dropbox/')
cor = readRDS('./metpred/dl/success_wgbs_0.2_10k/res/cor.rds')
str(cor)
samplecor = cor[[1]]
sitecor = cor[[2]]
meansamplecor = cor[[3]]
str(samplecor)
str(sitecor)
str(meansamplecor)
identical(names(samplecor), names(meansamplecor))
df = data.frame(correlation = c(sitecor, samplecor), type=c(rep('site',length(sitecor)), rep('sample',length(samplecor))))
library(ggplot2)
pdf(paste0('./metpred/dl/success_wgbs_0.2_10k/plot/sample_and_site_cor.pdf'),width=3,heigh=2)
ggplot(data=df,aes(x=type, y=correlation,fill=type)) + geom_violin(alpha=0.2,trim=T)  + 
  geom_boxplot(width=0.05,outlier.shape = NA) +
  theme_classic() +theme(legend.position = 'none', axis.text.x=element_text(color='black',size=11), axis.title = element_text(size=11)) + ylab('Spearman correlation')+
  xlab('')
dev.off()
