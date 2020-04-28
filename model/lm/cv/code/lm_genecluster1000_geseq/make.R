cvid <- as.numeric(commandArgs(trailingOnly = T))
load(paste0('/home-4/zji4@jhu.edu/scratch/metpred/model/lm/cv/res/lm_genecluster1000_geonly/',cvid,'.rda'))

seqclu <- readRDS('/home-4/zji4@jhu.edu/scratch/metpred/data/tcga/seq/final.rds')
seqclu <- seqclu[seqclu[,2] %in% row.names(am),]
tab <- table(seqclu[,1])
seqclu <- seqclu[seqclu[,1] %in% as.numeric(names(tab)[tab > 1]),]
neim <- t(sapply(unique(seqclu[,1]),function(i) colMeans(am[seqclu[seqclu[,1]==i,2],])))

predm <- neim[,trainid] %*% (trainx %*% chol2inv(chol(crossprod(trainx)))) %*% t(testx)

predm <- exp(predm)
predm <- predm/(1+predm)

predm <- predm[match(seqclu[,1],unique(seqclu[,1])),]
row.names(predm) <- seqclu[,2]

pred[row.names(predm),] <- (pred[row.names(predm),] + predm[row.names(predm),])/2
#rowcor <- sapply(1:nrow(testm),function(i) cor(testm[i,],pred[i,]))
#colcor <- sapply(1:ncol(testm),function(i) cor(testm[,i],pred[,i]))
rm('am')
save.image(file=paste0('/home-4/zji4@jhu.edu/scratch/metpred/model/lm/cv/res/lm_genecluster1000_geseq/',cvid,'.rda'))

