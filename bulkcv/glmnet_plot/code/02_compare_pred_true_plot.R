setwd('/home-4/whou10@jhu.edu/scratch/Wenpin/')
source('./resource/function.R')
load('./metpred/bulkcv/glmnet/arraypred.rda')
sitecor = sapply(1:nrow(pred),function(i){
  cor(pred[i,], true[i, ], method='spearman')
})
spcor = sapply(1:ncol(pred), function(i){
  cor(pred[,i], true[,i], method='spearman')
})
sitecor = unlist(sitecor)
spcor = unlist(spcor)
df = data.frame(correlation = c(sitecor, spcor), type = c(rep('samplewise', length(sitecor)), rep('sitewise',length(spcor))))
saveRDS(df,'./metpred/bulkcv/glmnet_plot/plot/compare_pred_true_pd.rds')
library(ggplot2)
pdf('./metpred/bulkcv/glmnet_plot/plot/compare_pred_true.pdf',width=3.5,height=2)
ggplot(data=df,aes(x=type, y=correlation,fill=type)) + geom_violin(alpha=0.2,trim=T)  + 
  geom_boxplot(width=0.1,outlier.shape = NA) +
  theme_classic() + xlab('Comparison between predicted and true')+theme(legend.position = 'none', axis.text.x=element_text(color='black',size=11), axis.title = element_text(size=11)) + ylab('Spearman correlation')+
  scale_fill_manual(values=c('blue','orange'))
dev.off()
