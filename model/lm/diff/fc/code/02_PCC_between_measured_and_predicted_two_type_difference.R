allf <- list.files('/home-4/whou10@jhu.edu/scratch/Wenpin/metpred/model/lm/diff/fc/tmpres/')
f <- allf[1]
cc <- sapply(allf, function(f){
  tmp <- readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/metpred/model/lm/diff/fc/tmpres/',f))
  return(tmp)
})

cc <- do.call(rbind, cc)

library(ggplot2)
library(reshape2)
library(RColorBrewer)
pd = melt(cc)
pdf('/home-4/whou10@jhu.edu/scratch/Wenpin/metpred/model/lm/diff/fc/plot/PCC_between_measured_and_predicted_two_type_difference.pdf',width=7,height=4)
ggplot(data=pd, aes(x=Var2, y=value, fill=Var2)) + geom_violin(alpha=0.1, scale='count') +
  geom_boxplot(width=0.05,alpha=0.5, outlier.shape = NA) +
  theme_classic() + xlab('') + ylab('Pearson correlation coefficient(PCC)') + ggtitle('PCC between measured and predicted two type difference') +
  theme(legend.position = 'none', axis.text.x=element_text(angle=45, hjust=1,  size=11, color='black'), axis.title.y = element_text(size=11, color='black')) +
  scale_fill_brewer(palette="Dark2") 
dev.off()
