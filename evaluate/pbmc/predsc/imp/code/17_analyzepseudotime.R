# library(BSgenome.Hsapiens.UCSC.hg38)
# k <- sapply(sample(1:nrow(me),40),function(i) {
#   as.character(getSeq(Hsapiens,sub('_.*','',rownames(me)[i]),start=as.numeric(sub('.*_','',rownames(me)[i])),width=2))
# })

suppressMessages(library(GenomicRanges))
me <- readRDS('/home/whou10/data/whou/metpred/evaluate/pbmc/predsc/imp/res/pseudotime/metest.rds')
ge <- readRDS('/home/whou10/data/whou/metpred/evaluate/pbmc/predsc/imp/res/pseudotime/getest.rds')
gtf <- data.table::fread('/home/whou10/scratch16/whou10/resource/gtf/grch38.gtf',data.table=F)
gtf <- gtf[gtf[,3]=='gene',]
gn <- sub('".*','',sub('.*gene_name "','',gtf[,9]))
gr <- GRanges(seqnames=gtf[,1],IRanges(start=gtf[,4],end=gtf[,5]),strand=gtf[,7])
names(gr) <- gn
gr <- gr[names(gr) %in% rownames(ge)]
pro <- promoters(gr,1000,0)
loc <- as.numeric(sub('.*_','',rownames(me)))
megr <- GRanges(seqnames=sub('_.*','',rownames(me)),IRanges(start=loc,end=loc))
names(megr) <- rownames(me)
o <- as.matrix(findOverlaps(pro,megr))
cor(ge[match(names(pro)[o[,1]],rownames(ge)),'fstat'],me[match(names(megr)[o[,2]],rownames(me)),'fstat'],method='spearman')

## =================================
gefit = readRDS('/home/whou10/data/whou/metpred/evaluate/pbmc/predsc/imp/res/pseudotime/gefit.rds')
mefit = readRDS('/home/whou10/data/whou/metpred/evaluate/pbmc/predsc/imp/res/pseudotime/mefit.rds')
me <- readRDS('/home/whou10/data/whou/metpred/evaluate/pbmc/predsc/imp/res/pseudotime/metest.rds')
mefit = mefit[rownames(me)[me[,1]<0.05],]


source('/home/whou10/scratch16/whou10/resource/myfunc/01_function.R')
## CpG in promoter region  <-> gene
## shorter distance, 

## find CpG in promoters
gtf <- data.table::fread('/home/whou10/scratch16/whou10/resource/gtf/grch38.gtf',data.table=F)
gtf <- gtf[gtf[,3]=='gene',]
gn <- sub('".*','',sub('.*gene_name "','',gtf[,9]))
library(GenomicRanges)
gr <- GRanges(seqnames=gtf[,1],IRanges(start=gtf[,4],end=gtf[,5]),strand=gtf[,7])
names(gr) <- gn
gr <- gr[names(gr) %in% rownames(gefit)]
pro <- promoters(gr,1000,0)
loc <- as.numeric(sub('.*_','',rownames(mefit)))
megr <- GRanges(seqnames=sub('_.*','',rownames(mefit)),IRanges(start=loc,end=loc))
names(megr) <- rownames(mefit)
o <- as.matrix(findOverlaps(pro,megr))

library(dtw)
dist <- sapply(unique(o[,1]),function(i) {
  print(i)
  mme <- colSums(mefit[names(megr)[o[o[,1]==i,2]],,drop=F])
  dtw_result <- dtw(gefit[names(pro)[i],],mme)
  alignment_distance <- dtw_result$distance
})
saveRDS(dist, '/home/whou10/data/whou/metpred/evaluate/pbmc/predsc/imp/res/pseudotime/mefitSum_gefit_dtw.rds')

set.seed(1)
dist2 <- sapply(sample(1:nrow(gefit), length(unique(o[,1]))),function(i) {
  print(i)
  mme <- colSums(mefit[sample(1:nrow(mefit), 10),, drop=F])
  dtw_result <- dtw(gefit[i,],mme)
  alignment_distance <- dtw_result$distance
})
saveRDS(dist2, '/home/whou10/data/whou/metpred/evaluate/pbmc/predsc/imp/res/pseudotime/mefitSum_gefit_dtw_random.rds')

data <- data.frame(
  value = c(dist, dist2),
  group = c(rep("Promoter", length(dist)), rep("Random", length(dist2)))
)

# Plot the overlaid histograms with transparency for each group
pdf('/home/whou10/data/whou/metpred/evaluate/pbmc/predsc/imp/plot/promoterDNAm_expr_dtw.pdf', 
    width = 3.1, height = 1.6)
ggplot(data, aes(x = group, y = value, fill = group)) +
  geom_boxplot(outliers = F, alpha = 0.7) +
  labs(title = "", y = "DTW dist (DNAm, expression)", x = "CpG region") +
  scale_fill_manual(values = c("chocolate", "darkgreen"))
dev.off()   

pdf('/home/whou10/data/whou/metpred/evaluate/pbmc/predsc/imp/plot/promoterDNAm_expr_dtw_hist.pdf', 
    width = 3.1, height = 1.6)
ggplot(data, aes(x = value, fill = group)) +
  geom_histogram(position = "identity", alpha = 0.5, bins = 30) +
  labs(title = "", x = "DTW dist (DNAm, expression)", y = "Count") +
  scale_fill_manual(values = c("skyblue", "lightcoral")) +
  xlim(0, 2e4) 
dev.off()

## =======================
## show an example gene 
## =======================
g = 'HLA-A' 
i = which(rownames(gefit)==g)
par(mfrow=c(1,2))
plot(1:ncol(gefit), gefit[i, ], col = 'red')
plot(colSums(mefit[o[o[,1]==i,2], ]))

pd = rbind(data.frame(value = gefit[g, ],
                      pseudotime = seq(1,ncol(gefit)),
                      type = 'expression'),
           data.frame(value = colSums(mefit[o[o[,1]==i,2], ]),
                pseudotime = seq(1,ncol(gefit)),
                type = 'promoter DNAm'))

pdf('/home/whou10/data/whou/metpred/evaluate/pbmc/predsc/imp/plot/promoterDNAm_expr_example.pdf', 
    width = 3.1, height = 1.6)
ggplot(data = pd, aes(x = pseudotime, y = value, color = type)) +
  geom_line() + 
  xlab('Pseudotime') + 
  ylab('Value(log2(TPM) or beta value)')+
  ggtitle(g)
dev.off()


library(pheatmap)
pdf('/home/whou10/data/whou/metpred/evaluate/pbmc/predsc/imp/plot/promoterDNAm_expr_example_hm.pdf', 
    width = 4.5, height = 2.5)
pheatmap(mefit[o[o[,1]==i,2], ], show_colnames = F, show_rownames = T,
         cluster_cols = F,
         scale = 'none', main = g)
dev.off()

## ==========================================
## show all promoter DNAM and gene expression
## ==========================================
ag = names(pro)[unique(o[,1])]
gefitsig = gefit[ag, ,drop=F]
str(gefitsig)
mefitsig <- t(sapply(unique(o[,1]), function(j){
  tmp = colSums(mefit[o[o[,1]==j,2], ,drop=F])
}))
str(mefitsig)

tmp = 
changepoint = sapply(1:nrow(gefitsig), function(i) {
  v <- scale(gefitsig)[i, seq(1, ncol(gefitsig) / 2)]
  which(v[-length(v)] * v[-1] < 0)[1]
})
summary(changepoint)
