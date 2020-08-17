source('/home-4/whou10@jhu.edu/scratch/Wenpin/metpred/software/predict.R')
source('/home-4/whou10@jhu.edu/scratch/Wenpin/metpred/software/trainmodel.R')
meth <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/metpred/data/combine/wgbs/sampleme.rds')
expr <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/metpred/data/combine/wgbs/filterge.rds')
genenametb <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/resource/hg38_geneid_genename_genelength.rds')
genenametb[,1] <- gsub('\\..*','', genenametb[,1])
rownames(expr) <- paste0(genenametb[match(rownames(expr), genenametb[,1]) , 2], ';', rownames(expr))
m <- trainmodel(expr, meth)
saveRDS(m, '/home-4/whou10@jhu.edu/scratch/Wenpin/metpred/model/wgbs/model.rds')


