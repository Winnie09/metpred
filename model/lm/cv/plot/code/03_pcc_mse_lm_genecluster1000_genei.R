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
  
allf = list.files('/home-4/zji4@jhu.edu/scratch/metpred/model/lm/cv/res/lm_genecluster1000_genei/')
pl_m <- pl_m2 <- pl <- pl2 <- list()
for (f in allf){
  print(f)
  load(paste0('/home-4/zji4@jhu.edu/scratch/metpred/model/lm/cv/res/lm_genecluster1000_genei/',f))
  model = model_perf(testm, pred)
  use_mean =  mean_perf(testm, trainy)
  pl[[f]] = melt(rbind(data.frame(model = model$crosssample, use_mean = use_mean$crosssample, measure='crosssample'),
  data.frame(model = model$crosssite, use_mean = use_mean$crosssite, measure='crosssite')))
  pl2[[f]] = melt(rbind(data.frame(model = model$crosssample_mse, use_mean = use_mean$crosssample_mse, measure='crosssample_mse'),
  data.frame(model = model$crosssite_mse, use_mean = use_mean$crosssite_mse, measure='crosssite_mse')))
  
  model = model_perf(testm[rownames(predm),], predm)
  pl_m[[f]] = melt(rbind(data.frame(model = model$crosssample, use_mean = use_mean$crosssample, measure='crosssample'),
    data.frame(model = model$crosssite, use_mean = use_mean$crosssite, measure='crosssite')))
  use_mean = mean_perf(neim[, colnames(testm)], neim[, colnames(trainy)])
  pl_m2[[f]] = melt(rbind(data.frame(model = model$crosssample_mse, use_mean = use_mean$crosssample_mse, measure='crosssample_mse'),
  data.frame(model = model$crosssite_mse, use_mean = use_mean$crosssite_mse, measure='crosssite_mse')))
  
}
  

pd = do.call(rbind, pl)
pd2 = do.call(rbind, pl2)
pd_nei = do.call(rbind, pl_m)
pd_nei2 = do.call(rbind, pl_m2)
saveRDS(list(pd=pd, pd2=pd2, pd_nei=pd_nei, pd_nei2=pd_nei2),'/home-4/whou10@jhu.edu/scratch/Wenpin/metpred/model/lm/plot/lm_genecluster1000_genei/pcc_mse_pd.rds')

p1 <- ggplot() + geom_boxplot(data=pd, aes(x=measure, y = value, fill=variable), alpha=0.2, outlier.size=0.1) + theme_classic() + xlab('') + ylab('Pearson correlation coefficient(PCC)') + ggtitle('average of locus-direct and neighboring-smoothing')
p2 <- ggplot() + geom_boxplot(data=pd2, aes(x=measure, y = value, fill=variable), alpha=0.2, outlier.size = 0.1) + theme_classic() + xlab('') + ylab('Mean Square Error(MSE)')

p3 <- ggplot() + geom_boxplot(data=pd_nei, aes(x=measure, y = value, fill=variable), alpha=0.2, outlier.size=0.1) + theme_classic() + xlab('') + ylab('Pearson correlation coefficient(PCC)') + ggtitle('Neighboring-smoothing only')
p4 <- ggplot() + geom_boxplot(data=pd_nei2, aes(x=measure, y = value, fill=variable), alpha=0.2, outlier.size = 0.1) + theme_classic() + xlab('') + ylab('Mean Square Error(MSE)')

ggsave('/home-4/whou10@jhu.edu/scratch/Wenpin/metpred/model/lm/plot/lm_genecluster1000_genei/pcc_mse.png',
       grid.arrange(p1,p2,p3,p4,nrow=2),
       width=8,height=6,dpi=200)

round(c(median(pd[pd$variable=='model' & pd$measure=='crosssample', 'value']),
median(pd[pd$variable=='use_mean' & pd$measure=='crosssample', 'value']),
median(pd[pd$variable=='model' & pd$measure=='crosssite', 'value']),
  median(pd[pd$variable=='use_mean' & pd$measure=='crosssite', 'value'])),2)

round(c(median(pd2[pd2$variable=='model' & pd2$measure=='crosssample_mse', 'value']),
median(pd2[pd2$variable=='use_mean' & pd2$measure=='crosssample_mse', 'value']),
median(pd2[pd2$variable=='model' & pd2$measure=='crosssite_mse', 'value']),
median(pd2[pd2$variable=='use_mean' & pd2$measure=='crosssite_mse', 'value'])),2)

