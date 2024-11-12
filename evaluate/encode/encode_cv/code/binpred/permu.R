ge <- readRDS('/home/whou10/data/whou/metpred/data/encode_wgbs/processed/GRCh38_normal/rna.rds')
n <- colnames(ge)
ge <- ge[,sample(1:ncol(ge))]
colnames(ge) <- n
me <- readRDS('/home/whou10/data/whou/metpred/data/encode_wgbs/processed/GRCh38_normal/wgbs_bin.rds')

# new program
source('/home/whou10/data/whou/metpred/software/trainpredict.R')
set.seed(12345)
cid <- as.numeric(cut(1:ncol(me),10))
names(cid) <- sample(colnames(me))

i <- as.numeric(commandArgs(trailingOnly = T))
testid <- which(cid==i)
trainid <- which(cid!=i)
pred <- trainpredict(trainexpr=ge[,trainid],testexpr=ge[,testid],trainmeth=me[,trainid])   
saveRDS(pred,file=paste0('/home/whou10/data/whou/metpred/evaluate/encode_cv/binpred/permu/',i,'.rds'))


# samplecv <- sapply(1:ncol(pred),function(i) cor(pred[,i],testme[,i]))
# sitecv <- sapply(1:nrow(pred),function(i) cor(pred[i,],testme[i,]))
# testsd <- apply(testme,1,sd)



