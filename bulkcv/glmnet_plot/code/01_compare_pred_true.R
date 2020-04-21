setwd('/home-4/whou10@jhu.edu/scratch/Wenpin/')
source('./resource/function.R')
load('./metpred/bulkcv/glmnet/arraypred.rda')
pred <- apply(pred,2,rank)
true <- apply(true,2,rank)
sitecor = corfunc(pred, true)
###
set.seed(12345)
id = sample( nrow(pred), 1e4)
pred = pred[id,]
true = true[id,]
rm = rm[id]
###
spcor = corfunc(t(pred),t(true))
spcor[lower.tri(spcor)] <- NA
sitecor[lower.tri(sitecor)] = NA
sitecor = as.vector(sitecor)
spcor = as.vector(spcor)


df = data.frame(correlation = c(sitecor, spcor), type = c(rep('samplewise', length(sitecor)), rep('sitewise',length(spcor))))
saveRDS(df,'./metpred/bulkcv/glmnet_plot/plot/compare_pred_true_pd.rds')
library(ggplot2)
pdf('./metpred/bulkcv/glmnet_plot/plot/compare_pred_true.pdf',width=3.5,height=2.5)
ggplot(data=df,aes(x=type, y=correlation,fill=type)) + geom_violin(alpha=0.2,trim=T)  + 
  geom_boxplot(width=0.1,outlier.shape = NA) +
  theme_classic() + xlab('Comparison between predicted and true')+theme(legend.position = 'none', axis.text.x=element_text(color='black',size=11), axis.title = element_text(size=11)) + ylab('Spearman correlation')+
  scale_fill_manual(values=c('blue','orange'))
dev.off()
