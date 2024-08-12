setwd('/home/whou10/data/whou/metpred/')
t = readRDS('./evaluate/eff/res/time.rds')
t = as.numeric(t[-1])
score = matrix(0,nrow=length(t),ncol=1)
rownames(score) = names(t)

max = max(t,na.rm=T)
min = min(t,na.rm=T)
score = 1- (t-min)/(max-min)
saveRDS(score,'./evaluate/eff/res/time_scaled.rds')
saveRDS(mean(score),'./evaluate/eff/res/time_finalscore.rds')

library(reshape2)
t = readRDS('./evaluate/eff/res/time.rds')
t = as.numeric(t[-1])
names(t) <- seq(1e2, 1e3+200, 100)
timedata = data.frame(sampleNumber = log10(as.numeric(names(t))), 
                      time = t,
                      stringsAsFactors = F)
v <- summary(lm(time~sampleNumber,data=timedata))$coefficients[2,1]  
saveRDS(v,'./evaluate/eff/res/scalability.rds')

tb = readRDS('./evaluate/eff/res/memory.rds')
tb = sub('K','',tb)
tb = as.numeric(tb[-1])
score = matrix(0,nrow=length(tb),ncol=1)
rownames(score) = names(tb)
max = max(tb,na.rm=T)
min = min(tb,na.rm=T)
score = 1- (tb-min)/(max-min)
saveRDS(score,'./evaluate/eff/res/memory_scaled.rds')
saveRDS(mean(score),'./evaluate/eff/res/memory_finalscore.rds')


