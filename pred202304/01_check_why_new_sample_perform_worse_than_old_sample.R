## working on ro rstudio_co 20230414
#new data
oge <- readRDS('/home/whou10/data/whou/metpred/data/encode_wgbs/processed/GRCh38_rna.rds') ## new
ome <- readRDS('/home/whou10/data/whou/metpred/data/encode_wgbs/processed/GRCh38_wgbs.rds')

#old data
ge <- readRDS('/home/whou10/data/whou/metpred/data/notused/encode_wgbs_version2020_notused/data/expr.rds') ## not used --> old data
me <- readRDS('/home/whou10/data/whou/metpred/data/notused/encode_wgbs_version2020_notused/data/wgbs.rds')

rownames(me)=sub('_',':',rownames(me))
colnames(ome) <- gsub(' ','_',sub('Homo sapiens ','',colnames(ome)))


intr <- intersect(rownames(ome),rownames(me))
intc <- intersect(colnames(ome),colnames(me))

omeint <- ome[intr,intc]
meint <- me[intr,intc]

v1 <- as.vector(omeint)
v2 <- as.vector(meint)
id <- which(!is.na(omeint))
cor(v1[id],v2[id])

me <- me[rowMeans(is.na(me))==0,]
me <- me[sample(1:nrow(me),1000000),]

# new program
source('/home/whou10/data/whou/metpred/software/trainpredict.R')
cid <- as.numeric(cut(1:ncol(me),3))
names(cid) <- sample(colnames(me))

i <- 3
testid <- names(which(cid==i))
trainid <- setdiff(colnames(me),testid)
pred <- trainpredict(trainexpr=ge[,trainid],testexpr=ge[,testid],trainmeth=me[,trainid])   
saveRDS(pred, '/home/whou10/data/whou/metpred/pred202304/pred.rds')

testme <- me[,testid]  
samplecv <- sapply(1:ncol(pred),function(i) cor(pred[,i],testme[,i]))
sitecv <- sapply(1:nrow(pred),function(i) cor(pred[i,],testme[i,]))
testsd <- apply(testme,1,sd)

summary(samplecv)

summary(sitecv)

summary(sitecv[testsd > 0.2])


# old program
source('/home/whou10/data/whou/metpred/software/pretraining/trainmodel.R')
source('/home/whou10/data/whou/metpred/software/pretraining/predict.R')

mod <- trainmodel(ge[,trainid],me[,trainid])
pred2 <- predict(ge[,testid],mod)

samplecv2 <- sapply(1:ncol(pred2),function(i) cor(pred2[,i],testme[,i]))
sitecv2 <- sapply(1:nrow(pred2),function(i) cor(pred2[i,],testme[i,]))

summary(samplecv)

summary(samplecv2)

summary(sitecv)

summary(sitecv2)

summary(sitecv[testsd > 0.2])

summary(sitecv2[testsd > 0.2])

## new code, old data, perf good
## new code, new data, not work
## use new code only
## new data, subset old sample, (increase sample, data processing is encode 4 might change the data values)
