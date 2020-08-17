source('/home-4/whou10@jhu.edu/scratch/metpred/software/predict.R')
source('/home-4/whou10@jhu.edu/scratch/metpred/software/trainmodel.R')
meth <- readRDS('/home-4/whou10@jhu.edu/scratch/metpred/data/combine/wgbs/sampleme.rds')
expr <- readRDS('/home-4/whou10@jhu.edu/scratch/metpred/data/combine/wgbs/filterge.rds')
ds <- sub(':.*','',colnames(meth))

for (testds in c('encode','blueprint','CEEHRC')) {
print(testds)
trainds <- setdiff(ds,testds)
m <- trainmodel(expr[,ds %in% trainds],meth[,ds %in% trainds])
pred <- predict(expr[,ds %in% testds],m)

scalematrix <- function(data) {
  cm <- rowMeans(data)
  csd <- sqrt((rowMeans(data*data) - cm^2) / (ncol(data) - 1) * ncol(data))
  (data - cm) / csd
}

corfunc <- function(m1,m2,type='concordant') {
  if (type=='concordant') {
    rowSums(scalematrix(m1) * scalematrix(m2))/(ncol(m1)-1)
  } else {
    scalematrix(t(m1)) %*% t(scalematrix(t(m2)))/(nrow(m1)-1)            
  }
}
cv <- corfunc(pred,meth[,ds %in% testds])
print(summary(cv))
cv <- corfunc(t(pred),t(meth[,ds %in% testds]))
print(summary(cv))
}
