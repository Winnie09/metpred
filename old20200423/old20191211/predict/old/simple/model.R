suppressMessages(library(GenomicRanges))
beta <- readRDS('/home-4/zji4@jhu.edu/scratch/metpred/data/mearray/beta.rds')

e <- readRDS('/home-4/zji4@jhu.edu/scratch/metpred/data/expr/mat.rds')
ecl <- readRDS('/home-4/zji4@jhu.edu/scratch/metpred/data/expr/clusterexpr.rds')
near <- readRDS('/home-4/zji4@jhu.edu/scratch/metpred/data/nearesttss/nearesttss.rds')
neargene <- e[near[,2],colnames(beta)]
neardist <- near[,3]

#kmer <- readRDS('/home-4/zji4@jhu.edu/scratch/metpred/data/kmer/flank3bp.rds')
bmean <- readRDS('/home-4/zji4@jhu.edu/scratch/metpred/data/mearray/beta_mean.rds')
bsd <- readRDS('/home-4/zji4@jhu.edu/scratch/metpred/data/mearray/beta_sd.rds')

library(xgboost)
set.seed(12345)
trainid <- sample(colnames(beta),110)
testid <- setdiff(colnames(beta),trainid)

mod <- NULL
for (sid in trainid) {
  print(sid)
  trainy <- as.vector(beta[,sid])
  trainx <- cbind(expr=neargene[,sid],neardist,bmean,bsd)
  if (is.null(mod)) {
    mod <- xgboost(trainx,trainy,objective='binary:logistic',nrounds=10,nthread=10)  
  } else {
    mod <- xgboost(trainx,trainy,objective='binary:logistic',nrounds=10,nthread=10,xgb_model=mod)
  }
}
saveRDS(mod,file='/home-4/zji4@jhu.edu/scratch/metpred/res/predict/model/mod.rds')
saveRDS(testid,file='/home-4/zji4@jhu.edu/scratch/metpred/res/predict/model/testid.rds')

