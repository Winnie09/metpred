load('/home-4/zji4@jhu.edu/scratch/metpred/bulkcv/pred.rda')
setwd('/home-4/whou10@jhu.edu/scratch/Wenpin/metpred/bulkcv')
sampcv <- sapply(1:ncol(pred),function(i) {
  cor(pred[,i],true[,i])
})
ctrlsampcv <- sapply(1:ncol(pred),function(i) {
  cor(rm,true[,i])
})

pdf('./acrosssite_correlation.pdf')
hist(sampcv, main='Within sample across sites Pearson correlation')
dev.off()

pdf('./acrosssite_correlation2.pdf')
plot(sampcv,ctrlsampcv,pch=20, main='Within sample across sites Pearson correlation',xlab='cor',ylab='frequency')
abline(0,1,col='red')
dev.off()


sitecv <- sapply(1:nrow(pred),function(i) {
  cor(pred[i,],true[i,])
})

pdf('./acrosssample_correlation.pdf')
hist(sitecv, main='Each site sample-wise Pearson correlation',xlab='cor',ylab='frequency')
abline(v=0,col='red')
dev.off()

pdf('./acrosssite_correlation_compare.pdf')
boxplot(data.frame(prediction=sampcv,meanprofile=ctrlsampcv))
dev.off()

