cvid <- as.numeric(commandArgs(trailingOnly = T))
load(paste0('/home-4/zji4@jhu.edu/scratch/metpred/model/lm/cv/res/lm_genecluster1000_geonly/',cvid,'.rda'))

nei <- readRDS('/home-4/zji4@jhu.edu/scratch/metpred/data/tcga/neighbor/10000.rds')
nei <- nei[sapply(nei,length) > 1]
nei <- nei[row.names(am)]
nei <- sapply(nei,function(i) intersect(i,row.names(am)))
nei <- nei[sapply(nei,length) > 1]
neim <- t(sapply(nei,function(i) colMeans(am[i,])))

predm <- neim[,trainid] %*% (trainx %*% chol2inv(chol(crossprod(trainx)))) %*% t(testx)

predm <- exp(predm)
predm <- predm/(1+predm)

pred[row.names(predm),] <- (pred[row.names(predm),] + predm[row.names(predm),])/2
#rowcor <- sapply(1:nrow(testm),function(i) cor(testm[i,],pred[i,]))
#colcor <- sapply(1:ncol(testm),function(i) cor(testm[,i],pred[,i]))
rm('am')
save.image(file=paste0('/home-4/zji4@jhu.edu/scratch/metpred/model/lm/cv/res/lm_genecluster1000_genei/',cvid,'.rda'))

