d <- readRDS('/home-4/zji4@jhu.edu/scratch/metpred/cellline/data/data/metclu/submet.rds')
dcl <- readRDS('/home-4/zji4@jhu.edu/scratch/metpred/cellline/data/data/metclu/clu1000.rds')
d <- t(sapply(1:max(dcl),function(i) colMeans(d[dcl==i,])))
ecl <- readRDS('/home-4/zji4@jhu.edu/scratch/metpred/cellline/data/data/database/clumean.rds')

set.seed(12345)
library(xgboost)
tecl <- t(ecl)
testidv <- c(rep(1:7,each=12),rep(8:10,each=11))
names(testidv) <- sample(colnames(d))
paralist <- expand.grid(eta=c(0.1,0.3),gamma=c(0,1))
id <- 1:nrow(d)
mod <- lapply(id,function(i) {
  perf <- sapply(1:nrow(paralist),function(pid) {
    mean(sapply(1:10,function(tid) {
      cvid <- names(testidv)[testidv==tid]
      trainid <- setdiff(names(testidv),cvid)
      mod <- xgboost(tecl[trainid,],d[i,trainid],object='reg:logistic',eta=unlist(paralist[pid,1]),gamma=unlist(paralist[pid,2]),nrounds=30,verbose=F)
      mean((predict(mod,tecl[cvid,])-d[i,cvid])^2)
    }))  
  })
  pid <- which.min(perf)
  xgboost(tecl,d[i,],object='reg:logistic',eta=unlist(paralist[pid,1]),gamma=unlist(paralist[pid,2]),nrounds=50,verbose=F)
})
saveRDS(mod,file=paste0('/home-4/zji4@jhu.edu/scratch/metpred/cellline/model/res/combinesite_logistic_tune/comb.rds'))


