expr = readRDS('/home/whou10/data/whou/metpred/software_test/rna_encode_sub_20.rds')
expr = expr[1:2e3,]
meth = readRDS('/home/whou10/data/whou/metpred/software_test/wgbs_encode_sub_20.rds')
meth[is.na(meth)] = 0
str(expr)
str(meth)

trainexpr <- expr[,1:11]
testexpr <- expr[,12:20]
trainmeth <- meth[,1:11]
testmeth <- meth[,12:20]

rm('expr')
rm('meth')

source('/home/whou10/data/whou/metpred/software/trainpredict.R')
pred <- trainpredict(trainexpr,testexpr,trainmeth) 

summary(sapply(1:ncol(pred),function(i) cor(pred[,i],testmeth[,i])))

summary(sapply(1:nrow(pred),function(i) cor(pred[i,],testmeth[i,])))

