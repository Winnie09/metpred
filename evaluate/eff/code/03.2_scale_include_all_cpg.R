rm(list=ls())
setwd('/home/whou10/data/whou/metpred/')
tb = readRDS('./evaluate/eff/res/time_include_allcpg.rds')
## =============================
## if >=1 mehtods
t = matrix(as.numeric(tb[,-1]),nrow=nrow(tb))
dimnames(t) = list(tb[,1],colnames(tb)[-1])
t.bak = t
t[is.na(t)] = 60*24*3 + 1
score = matrix(0,nrow=nrow(t),ncol=ncol(t))
dimnames(score) = dimnames(t)
for (i in 1:ncol(t)){
  max = max(t[,i],na.rm=T)
  min = min(t[,i],na.rm=T)
  score[,i] = 1- (t[,i]-min)/(max-min)
}
saveRDS(score,'./evaluate/eff/res/time_scaled_include_allcpg.rds')
saveRDS(rowMeans(score),'./evaluate/eff/res/time_finalscore_include_allcpg.rds')

library(reshape2)
t = t.bak
timedata <- melt(t)
colnames(timedata) <- c('method','sampleNumber','time')
timedata$sampleNumber <- log10(timedata$sampleNumber)
v <- sapply(unique(timedata$method),function(i){
  summary(lm(time~sampleNumber,data=timedata[timedata$method==i,]))$coefficients[2,1]  
})
names(v) = unique(timedata$method)
saveRDS(v,'./evaluate/eff/res/scalability_include_allcpg.rds')

tb = readRDS('./evaluate/eff/res/memory_include_allcpg.rds')
tb = sub('K','',tb)
t = matrix(as.numeric(tb[,-1]),nrow=nrow(tb))
dimnames(t) = list(tb[,1],colnames(tb)[-1])
t[is.na(t)] = max(t,na.rm = T)
score = matrix(0,nrow=nrow(t),ncol=ncol(t))
dimnames(score) = dimnames(t)
for (i in 1:ncol(t)){
  max = max(t[,i],na.rm=T)
  min = min(t[,i],na.rm=T)
  score[,i] = 1- (t[,i]-min)/(max-min)
}
saveRDS(score,'./evaluate/eff/res/memory_scaled_include_allcpg.rds')
saveRDS(rowMeans(score),'./evaluate/eff/res/memory_finalscore_include_allcpg.rds')

