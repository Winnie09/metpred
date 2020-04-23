load('/home-4/zji4@jhu.edu/scratch/metpred/bulkcv/pred.rda')
sampcv <- sapply(1:ncol(pred),function(i) {
  cor(pred[,i],true[,i])
})
pdf('/home-4/zji4@jhu.edu/scratch/metpred/bulkcv/acrosssite_correlation.pdf')
hist(sampcv)
dev.off()
sitecv <- sapply(1:nrow(pred),function(i) {
  cor(pred[i,],true[i,])
})
pdf('/home-4/zji4@jhu.edu/scratch/metpred/bulkcv/acrosssample_correlation.pdf')
hist(sitecv)
dev.off()

ctrlsampcv <- sapply(1:ncol(pred),function(i) {
  cor(rm,true[,i])
})

pdf('/home-4/zji4@jhu.edu/scratch/metpred/bulkcv/acrosssite_correlation_compare.pdf')
boxplot(data.frame(prediction=sampcv,meanprofile=ctrlsampcv))
dev.off()

