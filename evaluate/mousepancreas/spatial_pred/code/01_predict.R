## sd of 
expr <- readRDS('/home/whou10/data/whou/metpred/data/mousepancreas/spatial/normexpr.rds')
dim(expr)
alls <- sub(':.*', '', colnames(expr))
s = alls[1]
source('/home/whou10/scratch16/whou10/resource/startup.R')

rs <- sapply(unique(alls), function(s){
  e.s <- expr[, alls == s]
  rs = rowsds(e.s)
})
str(rs)  
library(pheatmap)
# pheatmap(rs,scale = 'none', show_rownames = F, show_colnames = T)

## input: training expression: 
## all others cell type specific gene expression, alll other cell type specific wgbs as train dnam
## one sample's 9k genes 

## leave out 1k genes: gene.sel
prome <-  readRDS('/home/whou10/data/whou/FeinbergLab/mousePancreas/wgbs/promoter_methylation/promoter_tssup1kb_averaged_dnam.rds')
int = intersect(rownames(prome), rownames(rs))
rs = rs[int, ]
prome = prome[int, ]
v = rs[,1]
names(v) = rownames(rs)
gene.sel <- names(head(sort(v, decreasing = F), 1e3))


## select test sample 
test.s <- unique(alls)[1]
train.s <- setdiff(colnames(rs), test.s)

## test samples' gs: prome
prome <- prome[which(rownames(prome) %in% gene.sel), ]
str(prome)
prome.cn <- colnames(prome)
prome.cn <- sub('_unrecovered', '', sub('_acinar', '', sub('_duct', '', sub('_ADM', '', prome.cn))))
prome <- prome[, prome.cn %in% test.s]

## load cell type specific expression and wgbs data
r <- readRDS('/home/whou10/data/whou/metpred/data/mousepancreas/spatial/cse.rds')
w <- readRDS('/home/whou10/data/whou/metpred/data/mousepancreas/wgbs/bs.rds')

## select the CpGs in 1k genes' promoters
suppressMessages(library(bsseq))
library(data.table)
# load('/home/whou10/data/whou/FeinbergLab/mousePancreas/wgbs/data/processed/adm_bsseq_smoothed_cov_filtered.rda')
#m <- assays(bsseq_smoothed_sub)@listData$M
#cov <- assays(bsseq_smoothed_sub)@listData$Cov
gtf <- fread('/home/whou10/scratch16/whou10/resource/gtf/grcm38.gtf',data.table=F)
gtf <- gtf[gtf[,3]=='gene',]
gene <- GRanges(seqnames=gtf[,1],IRanges(start=gtf[,4],end=gtf[,5]),strand=gtf[,7])
names(gene) <- sub('".*','',sub('.*gene_name "','',gtf[,9]))
gene <- gene[gene.sel]
pro <- promoters(gene,upstream=1000,downstream=0)
wgbs.cpg = GRanges(seqnames = sub(':.*', '', rownames(w)), 
                   IRanges(start = sub('.*:', '', rownames(w))))
# o <- as.matrix(findOverlaps(bsseq_smoothed_sub@rowRanges,pro))
o <- as.matrix(findOverlaps(wgbs.cpg,pro))
str(o)
cpg.id <- rownames(w)[o[,1]]
str(cpg.id)
length(unique(o[,2]))


########## train model and predict
int <- intersect(colnames(r),colnames(w))
samp <- sub(':.*','',int)
# > str(samp)
#  chr [1:41] "A29_M_P_K4_D2" "A29_M_P_K4_D2" "A34_M_P_K4_D2" "A34_M_P_K4_D2" ...
trainid <- int[which(samp != test.s)]
testid <- setdiff(int, trainid)
library(Matrix)
# e <- readRDS('/home/whou10/data/whou/metpred/data/mousepancreas/spatial/normexpr.rds')
se <- sub(':.*','',colnames(expr))
expr <- expr[, se == test.s]
str(expr)
source('/home/whou10/data/whou/metpred/software/trainpredict.R')
source('/home/whou10/scratch16/whou10/resource/startup.R')
pred <- trainpredict(trainexpr=r[!rownames(r) %in% gene.sel,trainid],
                     testexpr=expr[!rownames(expr) %in% gene.sel,],
                     trainmeth=w[cpg.id, trainid])
setwd('/home/whou10/data/whou/metpred/evaluate/mousepancreas/')
saveRDS(pred, paste0('spatial_pred/res/pred_', test.s, '.rds'))

