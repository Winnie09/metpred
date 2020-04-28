model_perf <- function(testm, pred){
  rowcor <- sapply(1:nrow(testm),function(i) cor(testm[i,],pred[i,]))
  colcor <- sapply(1:ncol(testm),function(i) cor(testm[,i],pred[,i]))
  se <- (testm-pred)^2
  rowmse <- rowMeans(se)
  colmse <- colMeans(se)
  model = list(crosssample = rowcor, crosssite = colcor, crosssample_mse = rowmse, crosssite_mse = colmse)
  return(model)
}
  
mean_perf <- function(testm, trainy){
  meany = rowMeans(trainy)
  rowcor <- 0
  colcor <- sapply(1:ncol(testm),function(i) cor(testm[,i], meany))
  se <- (testm-meany)^2
  rowmse <- rowMeans(se)
  colmse <- colMeans(se)
  use_mean = list(crosssample = rowcor, crosssite = colcor, crosssample_mse = rowmse, crosssite_mse = colmse)
  return(use_mean)
}

library(ggplot2)
library(gridExtra)
library(reshape2)
  
allf = list.files('/home-4/zji4@jhu.edu/scratch/metpred/model/lm/cv/res/lm_genecluster1000_geseq/')
pl <- pl2 <- list()
for (f in allf){
  load(paste0('/home-4/zji4@jhu.edu/scratch/metpred/model/lm/cv/res/lm_genecluster1000_geseq/',f))
  model = model_perf(testm, pred)
  use_mean =  mean_perf(testm, trainy)
  pl[[f]] = melt(rbind(data.frame(model = model$crosssample, use_mean = use_mean$crosssample, measure='crosssample'),
  data.frame(model = model$crosssite, use_mean = use_mean$crosssite, measure='crosssite')))
  pl2[[f]] = melt(rbind(data.frame(model = model$crosssample_mse, use_mean = use_mean$crosssample_mse, measure='crosssample_mse'),
  data.frame(model = model$crosssite_mse, use_mean = use_mean$crosssite_mse, measure='crosssite_mse')))
}
  
pd = do.call(rbind, pl)
pd2 = do.call(rbind, pl2)

p1 <- ggplot() + geom_boxplot(data=pd, aes(x=measure, y = value, fill=variable), alpha=0.2, outlier.size=0.1) + theme_classic() + xlab('') + ylab('Pearson correlation coefficient(PCC)')
p2 <- ggplot() + geom_boxplot(data=pd2, aes(x=measure, y = value, fill=variable), alpha=0.2, outlier.size = 0.1) + theme_classic() + xlab('') + ylab('Mean Square Error(MSE)')

saveRDS(list(pd=pd, pd2=pd2),'/home-4/whou10@jhu.edu/scratch/Wenpin/metpred/model/lm/cv/plot/lm_genecluster1000_geseq/pcc_mse_pd.rds')

ggsave('/home-4/whou10@jhu.edu/scratch/Wenpin/metpred/model/lm/cv/plot/lm_genecluster1000_geseq/pcc_mse.png',
       grid.arrange(p1,p2,nrow=1),
       width=8,height=3,dpi=200)

a <- round(c(median(pd[pd$variable=='model' & pd$measure=='crosssample', 'value']),
median(pd[pd$variable=='use_mean' & pd$measure=='crosssample', 'value']),
median(pd[pd$variable=='model' & pd$measure=='crosssite', 'value']),
median(pd[pd$variable=='use_mean' & pd$measure=='crosssite', 'value'])),2)
print(a)

b <- round(c(median(pd2[pd2$variable=='model' & pd2$measure=='crosssample_mse', 'value']),
median(pd2[pd2$variable=='use_mean' & pd2$measure=='crosssample_mse', 'value']),
median(pd2[pd2$variable=='model' & pd2$measure=='crosssite_mse', 'value']),
median(pd2[pd2$variable=='use_mean' & pd2$measure=='crosssite_mse', 'value'])),2)
print(b)

