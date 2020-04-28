ap <- readRDS('/home-4/zji4@jhu.edu/scratch/metpred/data/tcga/proc/combine/project.rds')
am <- readRDS('/home-4/zji4@jhu.edu/scratch/metpred/data/tcga/proc/combine/completehighsme.rds')
ae <- readRDS('/home-4/zji4@jhu.edu/scratch/metpred/data/tcga/proc/combine/ge.rds')

spid <- split(sample(names(ap)),ceiling(seq_along(names(ap))/(length(ap)/10)))

cvid <- as.numeric(commandArgs(trailingOnly = T))
testid <- spid[[cvid]]
trainid <- setdiff(names(ap),testid)
e <- ae[,trainid]
e <- e[rowMeans(e >= 5) >= 0.01,]
m <- rowMeans(e)
s <- apply(e,1,sd)
cv <- s/m
id <- which(cv >= 0.1)
e <- e[id,]
m <- m[id]
s <- s[id]
stde <- (e-m)/s
clun <- 1000
set.seed(12345)
clu <- kmeans(stde,clun,iter.max = 10000)$cluster
clue <- sapply(1:clun,function(i) colMeans(e[names(clu)[clu==i],,drop=F]))
testclue <- sapply(1:clun,function(i) colMeans(ae[names(clu)[clu==i],testid,drop=F]))
trainx <- cbind(1,clue)
testx <- cbind(1,testclue)
trainy <- am[,trainid]
trainy[trainy==0] <- min(trainy[trainy>0])
trainy[trainy==1] <- max(trainy[trainy<1])
trainy <- log(trainy/(1-trainy))

# glmnet
# library(parallel)
# library(glmnet)
# pred <- mclapply(1:nrow(trainy), function(i) {
#   cvfit <- cv.glmnet(trainx,trainy[i,],family='gaussian')
#   predict(cvfit,newx=testx,s='lambda.min',type='response')[,1]
# },mc.cores=detectCores()-1)
# pred <- do.call(rbind,pred)

# lm
pred <- trainy %*% (trainx %*% chol2inv(chol(crossprod(trainx)))) %*% t(testx)


pred <- exp(pred)
pred <- pred/(1+pred)
testm <- am[,testid]
#rowcor <- sapply(1:nrow(testm),function(i) cor(testm[i,],pred[i,]))
#colcor <- sapply(1:ncol(testm),function(i) cor(testm[,i],pred[,i]))

save.image(file=paste0('/home-4/zji4@jhu.edu/scratch/metpred/model/lm/cv/res/',cvid,'.rda'))
