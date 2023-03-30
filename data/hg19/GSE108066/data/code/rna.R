library(data.table)
d1 <- fread('/home-4/zji4@jhu.edu/scratch/metpred/data/GSE108066/data/raw/rna/GSE107638_COUNTS_Exons_NeuN.txt.gz',data.table=F)
row.names(d1) <- d1[,1]
d2 <- fread('/home-4/zji4@jhu.edu/scratch/metpred/data/GSE108066/data/raw/rna/GSE107638_COUNTS_Exons_OLIG2.txt.gz',data.table=F)
row.names(d2) <- d2[,1]
d <- cbind(d1[,-1],d2[,-1])
af <- list.files('/home-4/zji4@jhu.edu/scratch/metpred/data/GSE108066/data/raw/wgbs')
af <- sub('GSM[0-9]*_','',sub('_CpG_WGBS.txt.gz','',af))
cn <- colnames(d)
cn <- sub('^X','',sub('_Schizo','',sub('_Control','',cn)))
colnames(d) <- cn
d <- d[,colnames(d) %in% af]
colnames(d) <- paste0('GSE108066_',colnames(d))
gl <- readRDS('/home-4/zji4@jhu.edu/scratch/resource/gl/res/hg19.rds')
load('/home-4/zji4@jhu.edu/scratch/resource/gn/res/hg19.rda')
d <- d[rownames(d) %in% genename,]
gl <- gl[geneid[match(rownames(d),genename)]]
d <- d/gl*10^3
d <- t(t(d)/colSums(d)*10^6)
d <- log2(d + 1)
saveRDS(d,file='/home-4/zji4@jhu.edu/scratch/metpred/data/GSE108066/data/proc/rna.rds')
