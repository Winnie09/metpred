suppressMessages(library(GenomicRanges))
load('/home-4/zji4@jhu.edu/scratch/encode_compiled/May19/metadata/DNase/processed/processed.rda')
af <- sub('.rds','',list.files('/home-4/zji4@jhu.edu/scratch/scate/data/hg19/res/DNasenormfeat/200bp'))
hg19_experiment <- hg19_experiment[hg19_experiment[,'Accession'] %in% af,]
dne <- gsub(' ','_',hg19_experiment[,'Biosample summary'])
beta <- readRDS('/home-4/zji4@jhu.edu/scratch/metpred/data/mearray/beta.rds')
#beta <- matrix(as.numeric(beta > 0.5),nrow=nrow(beta),dimnames = dimnames(beta))
d <- readRDS('/home-4/zji4@jhu.edu/scratch/encode_compiled/Nov19/mearray/data/mat/gr.rds')
gr <- readRDS('/home-4/zji4@jhu.edu/scratch/scate/data/hg19/res/genomegr/200bp/gr.rds')
o <- as.matrix(findOverlaps(d,gr))

amat <- sapply(intersect(dne,colnames(beta)),function(i) {
  tar <- hg19_experiment[dne==i,'Accession']
  readRDS(paste0('/home-4/zji4@jhu.edu/scratch/scate/data/hg19/res/DNasenormfeat/200bp/',tar,'.rds'))
})
library(preprocessCore)
dn <- dimnames(amat)
amat <- normalize.quantiles(amat)
dimnames(amat) <- dn
amat <- amat[o[,2],]

bmat <- beta[o[,1],intersect(dne,colnames(beta))]

near <- readRDS('/home-4/zji4@jhu.edu/scratch/metpred/data/nearesttss/nearesttss.rds')
e <- readRDS('/home-4/zji4@jhu.edu/scratch/metpred/data/expr/mat.rds')

neargene <- e[near[o[,1],2],colnames(amat)]
neardist <- near[o[,1],3]

library(xgboost)
set.seed(12345)
trainid <- sample(colnames(bmat),30)
testid <- setdiff(colnames(bmat),trainid)
trainy <- as.vector(bmat[,trainid])
testy <- as.vector(bmat[,testid])
trainx <- cbind(as.vector(amat[,trainid]),as.vector(neargene[,trainid]),neardist)
testx <- cbind(as.vector(amat[,testid]),as.vector(neargene[,testid]),neardist)
mod <- xgboost(trainx,trainy,objective='binary:logistic',nrounds=50)
pred <- predict(mod,testx)
#pred <- as.numeric(pred > 0.5)
#library(ROCR)
#perf <- prediction(pred,testy)
save.image('/home-4/zji4@jhu.edu/scratch/metpred/code/predict/res.rda')

#auc.perf = performance(perf, measure = "auc")
#pdf('/home-4/zji4@jhu.edu/scratch/metpred/code/predict/perf.pdf')
#plot(performance(perf, measure = "tpr", x.measure = "fpr"),main=round(auc.perf@y.values[[1]],4))
#dev.off()

#200bp 0.8070766

testy <- matrix(testy,ncol=length(testid))
pred <- matrix(pred,ncol=length(testid))

sampcv <- sapply(1:ncol(testy),function(s) {
  cor(testy[,s],pred[,s])
})

sitecv <- sapply(1:nrow(testy),function(s) {
  cor(testy[s,],pred[s,])
})

pdf('/home-4/zji4@jhu.edu/scratch/metpred/code/predict/sampcv.pdf')
hist(sampcv)
dev.off()

pdf('/home-4/zji4@jhu.edu/scratch/metpred/code/predict/sitecv.pdf')
hist(sitecv)
dev.off()
