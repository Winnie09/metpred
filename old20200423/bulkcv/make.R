# read in dna methylation data
d <- readRDS('/home-4/zji4@jhu.edu/scratch/metpred/data/data/mearray/beta.rds')
library(matrixStats)
# 20 as testing and remaining as training
testid <- sample(colnames(d),20)
trainid <- setdiff(colnames(d),testid)
true <- d[,testid]
d <- d[,trainid]
# only retain > 0.2 sd CpG sites
id <- which(rowSds(d) > 0.2)
true <- true[id,]
d <- d[id,]

m <- readRDS('/home-4/zji4@jhu.edu/scratch/metpred/data/data/expr/mat.rds')
m <- m[,trainid]
m <- m[rowSums(m >= 1) >= 1,]
library(preprocessCore)
dn <- dimnames(m)
m <- normalize.quantiles(m)
dimnames(m) <- dn
clu <- kmeans(t(apply(m,1,scale)),nrow(m)/100,iter.max=10000)$cluster
b <- m
ecl <- t(sapply(1:max(clu),function(i) {
  colMeans(m[names(clu)[clu==i],])
}))

# cluster and average CpG sites
scaled <- t(apply(d,1,scale))
dclu <- kmeans(scaled,1000,iter.max = 10000)$cluster
dmean <- t(sapply(1:max(dclu),function(i) {
  colMeans(d[dclu==i,])
}))

m <- readRDS('/home-4/zji4@jhu.edu/scratch/metpred/data/data/expr/mat.rds')
m <- m[,testid]
intg <- intersect(row.names(m),row.names(b))
m <- m[intg,]
b <- b[intg,]
bm <- rowMeans(apply(b,2,sort))
rn <- row.names(m)
m <- apply(m,2,function(i) bm[rank(i)])
row.names(m) <- rn
clu <- clu[names(clu) %in% row.names(m)]
scecl <- t(sapply(1:max(clu),function(i) {
  colMeans(m[names(clu)[clu==i],])
}))


library(glmnet)
tecl <- t(ecl)
colnames(tecl) <- 1:ncol(tecl)
tscecl <- t(scecl)
colnames(tscecl) <- 1:ncol(tscecl)

library(parallel)
pred <- mclapply(1:nrow(d), function(i) {
  # corv <- abs(apply(tecl,2,cor,d[i,]))
  # id <- which(corv >= sort(corv,decreasing = T)[10])
  # df <- data.frame(y=d[i,],tecl[,id])
  # fit <- glm(y~.,data=df,family='binomial',weights=rep(1,ncol(d)))
  # predict(fit,data.frame(tscecl[,id]),type='response')
  cvfit <- cv.glmnet(tecl,cbind(1-d[i,],d[i,]),family='binomial')
  predict(cvfit,newx=tscecl,s='lambda.min',type='response')[,1]
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
  cor(finalpred[,i],true[,i])
})

sitecv <- sapply(1:nrow(true),function(i) {
  cor(finalpred[i,],true[i,])
})
rm=rowMeans(d)
save(list=c('true','pred','rm'),file='/home-4/zji4@jhu.edu/scratch/metpred/bulkcv/pred.rda')

