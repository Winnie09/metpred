library(ggplot2)
library(gridExtra)

rdir <- '/home-4/whou10@jhu.edu/scratch/Wenpin/metpred/res/res/'

for (disease in c('cancer','normal')) {
  pred <- readRDS(paste0(rdir, 'predicted_DNAm_on_testset_', disease,'.rds'))
  cv <- readRDS(paste0(rdir, 'sampcvall_', disease, '.rds'))
  sdv <- apply(pred[names(cv),],1,sd)
  cc <- cut(sdv,seq(0,1,0.05))
  pdf(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/metpred/res/plot/predictedvaluesd_cor/',disease,'.pdf'),width=6,height=6)
  print(ggplot(data.frame(cor=cc,cv=cv),aes(y=cv,x=cor,fill=cor)) + geom_violin(alpha=0.3) + theme_classic() + scale_fill_manual(values=rainbow(length(unique(cc)))) + xlab('Predicted value sd') + ylab('Correlation') + theme(legend.position = 'none', axis.text.x = element_text(angle = 45, hjust = 1)))
  dev.off()
}


