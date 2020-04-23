trainx <- t(readRDS('/home-4/zji4@jhu.edu/scratch/metpred/bulkcv//data/trainx.rds'))
trainy <- readRDS('/home-4/zji4@jhu.edu/scratch/metpred/bulkcv//data/trainy.rds')
testx <- t(readRDS('/home-4/zji4@jhu.edu/scratch/metpred/bulkcv//data/testx.rds'))
testy <- readRDS('/home-4/zji4@jhu.edu/scratch/metpred/bulkcv//data/testy.rds')

library(parallel)
library(glmnet)
pred <- mclapply(1:nrow(trainy), function(i) {
  cvfit <- cv.glmnet(trainx,cbind(1-trainy[i,],trainy[i,]),family='binomial')
  predict(cvfit,newx=testx,s='lambda.min',type='response')[,1]
},mc.cores=detectCores()-1)

pred <- do.call(rbind,pred)

# meanpred <- mclapply(1:nrow(dmean), function(i) {
#   corv <- abs(apply(tecl,2,cor,dmean[i,]))
#   id <- which(corv >= sort(corv,decreasing = T)[10])
#   df <- data.frame(y=dmean[i,],tecl[,id])
#   fit <- glm(y~.,data=df,family='binomial',weights=rep(1,ncol(d)))
#   predict(fit,data.frame(tscecl[,id]),type='response')
# },mc.cores=detectCores()-1)
# meanpred <- do.call(rbind,meanpred)
# 
# finalpred <- (pred+meanpred[dclu,])/2
finalpred <- pred
# calculate correlation across sites or samples
samplecv <- sapply(1:ncol(pred),function(i) {
  cor(finalpred[,i],testy[,i])
})

sitecv <- sapply(1:nrow(pred),function(i) {
  cor(finalpred[i,],testy[i,])
})
#rm=rowMeans(d)
save.image(file='/home-4/zji4@jhu.edu/scratch/metpred/bulkcv/pred.rda')

