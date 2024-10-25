predge <- readRDS('/home/whou10/data/whou/metpred/data/gtex_wgbs/processed/rna.rds')
ge <- readRDS('/home/whou10/data/whou/metpred/data/encode_wgbs/processed/GRCh38_normal/rna.rds')
rownames(predge) <- sub(':.*','',rownames(predge))

int <- intersect(rownames(ge),rownames(predge))
ge <- ge[int,]
predge <- predge[int,]

me <- readRDS('/home/whou10/data/whou/metpred/data/encode_wgbs/processed/GRCh38_normal/wgbs_bin.rds')

source('/home/whou10/data/whou/metpred/software/trainpredict.R')
pred <- trainpredict(trainexpr=ge,testexpr=predge,trainmeth=me)
saveRDS(pred,file='/home/whou10/data/whou/metpred/evaluate/encode_gtex/binpred/ramp.rds')


#c1 <- sapply(1:nrow(predme),function(i) cor(predme[i,],pred[i,]))
#c2 <- sapply(1:ncol(predme),function(i) cor(predme[,i],pred[,i]))

#sd <- apply(predme,1,sd)

