setwd('/home/whou10/data/whou/metpred/evaluate/mousepancreas/')
dir.r <- 'spatialpred/res/'
source('/home/whou10/scratch16/whou10/resource/startup.R')

expr <- readRDS('/home/whou10/data/whou/metpred/data/mousepancreas/spatial/normexpr.rds')
alls <- sub(':.*', '', colnames(expr))
sid <- as.numeric(commandArgs(trailingOnly = T)[[1]])
print(sid)

s = unique(alls)[sid]
rs <- sapply(unique(alls), function(s){
  e.s <- expr[, alls == s]
  rs = rowsds(e.s)
})
# library(pheatmap)
# pheatmap(rs,scale = 'none', show_rownames = F, show_colnames = T)

## leave out 1k genes: gene.sel
prome <-  readRDS('/home/whou10/data/whou/FeinbergLab/mousePancreas/wgbs/promoter_methylation/promoter_tssup1kb_averaged_dnam.rds')
int = intersect(rownames(prome), rownames(rs))
rs = rs[int, ]
prome = prome[int, ]
v = rs[,s]
names(v) = rownames(rs)
gene.sel <- names(head(sort(v, decreasing = T), 1e3))

## select test sample 
test.s <- s
train.s <- setdiff(colnames(rs), test.s)

## load cell type specific expression and wgbs data
r <- readRDS('/home/whou10/data/whou/metpred/data/mousepancreas/spatial/cse.rds')
w <- readRDS('/home/whou10/data/whou/metpred/data/mousepancreas/wgbs/bs.rds')

## select the CpGs in 1k genes' promoters
suppressMessages(library(bsseq))
library(data.table)
gtf <- fread('/home/whou10/scratch16/whou10/resource/gtf/grcm38.gtf',data.table=F)
gtf <- gtf[gtf[,3]=='gene',]
gene <- GRanges(seqnames=gtf[,1],IRanges(start=gtf[,4],end=gtf[,5]),strand=gtf[,7])
names(gene) <- sub('".*','',sub('.*gene_name "','',gtf[,9]))
gene <- gene[gene.sel]
pro <- promoters(gene,upstream=1000,downstream=0)
wgbs.cpg = GRanges(seqnames = sub(':.*', '', rownames(w)), 
                   IRanges(start = sub('.*:', '', rownames(w))))
o <- as.matrix(findOverlaps(wgbs.cpg,pro))
cpg.id <- rownames(w)[o[,1]]
saveRDS(data.frame(cpg = rownames(w)[o[,1]], gene = names(pro)[o[,2]]), 
        paste0('spatialpred/genecpg/', s, '.rds'))

