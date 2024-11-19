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
    width = 2.2, height = 2.2)
ggplot(data, aes(x = group, y = value, fill = group)) +
  geom_boxplot(outliers = F, alpha = 0.7) +
  labs(title = "", y = "DTW dist (DNAm, expression)", x = "CpG region") +
  scale_fill_manual(values = c("chocolate", "darkgreen")) +
  theme(legend.position = 'bottom')
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
names(pro)[unique(o[,1])]
g = 'GADD45A' #'CD27' #'CD45RA' #'CCR7' # 'IL2RA' #'FOXP3'
# level1:'GADD45A' , 'FCER1G'
# level2: 'NRAS' 'FCRL3' 'SPI1' # 'TGFBR2'#  # 'IKZF3' #'STAT4' #'IKZF1' #'CD244' 
i = which(names(pro)==g)
o[o[,1]==i,2]
par(mfrow=c(1,2))
plot(1:ncol(gefit), gefit[i, ], col = 'red', ylim=c(0,5))
plot(colSums(mefit[o[o[,1]==i,2], ,drop=F]), ylim=c(-0.1,12))


pd = rbind(data.frame(value = scale(gefit[g, ]),
                      pseudotime = seq(1,ncol(gefit)),
                      type = 'expression'),
           data.frame(value = scale(colSums(mefit[o[o[,1]==i,2], ,drop=F])),
                      pseudotime = seq(1,ncol(gefit)),
                      type = 'promoter DNAm'))
saveRDS(pd, '/home/whou10/data/whou/metpred/evaluate/pbmc/predsc/imp/plot/promoterDNAm_expr_example.rds')
pdf('/home/whou10/data/whou/metpred/evaluate/pbmc/predsc/imp/plot/promoterDNAm_expr_example.pdf', 
    width = 3.1, height = 1.6)
ggplot(data = pd, aes(x = pseudotime, y = value, color = type)) +
  geom_line() + 
  xlab('Pseudotime') + 
  ylab('Scaled log2(TPM) or beta value')+
  ggtitle(g)
dev.off()


library(pheatmap)
pdf('/home/whou10/data/whou/metpred/evaluate/pbmc/predsc/imp/plot/promoterDNAm_expr_example_hm.pdf', 
    width = 4, height = 2.5)
pheatmap(mefit[o[o[,1]==i,2], ], show_colnames = F, show_rownames = T,
         cluster_cols = F,
         scale = 'none', main = g)
dev.off()

## ==========================================
## show all promoter DNAM and gene expression
## ==========================================
range(unique(o[,1]))
length(pro)
ag = names(pro)[unique(o[,1])]
gefitsig = gefit[ag, ,drop=F]
str(gefitsig)
mefitsig <- t(sapply(unique(o[,1]), function(j){
  tmp = colSums(mefit[o[o[,1]==j,2], ,drop=F])
}))
str(mefitsig)
rownames(mefitsig) = ag


changepoint = sapply(1:nrow(gefitsig), function(i) {
  v <- scale(gefitsig[i, ])
  # res = which(v[-length(v)] * v[-1] < 0)[1]
  # print(res)
  # res
  
  # Calculate the differences
  diff_vec <- diff(v)
  
  # Identify the first stationary point
  # Threshold is set to a small value (e.g., 1e-6) to detect near-zero changes
  threshold <- 1e-3
  stationary_index <- which(abs(diff_vec) < threshold)[1] + 1  # +1 to adjust for diff() reducing length
  
  # Output the first stationary point
  if (!is.na(stationary_index)) {
    # cat("The first stationary point is at index:", stationary_index, "\n")
    # cat("Value at the stationary point:", vec[stationary_index], "\n")
    return(stationary_index)
  } else {
    # cat("No stationary point found within the threshold.\n")
    return(ncol(gefitsig))
  }
}, USE.NAMES = T)
str(changepoint)
summary(changepoint)
names(changepoint) = rownames(gefitsig)

gefitsig = gefitsig[order(changepoint), ]
mefitsig = mefitsig[rownames(gefitsig),  ]
pheatmap_result <- pheatmap(gefitsig, cluster_rows = T, cluster_cols = F, show_rownames = T, show_colnames = F, scale = 'row')
row_clustering <- pheatmap_result$tree_row$order
pdf('/home/whou10/data/whou/metpred/evaluate/pbmc/predsc/imp/plot/gefit_hm.pdf', 
    width = 12, height = 100)
pheatmap_result
dev.off()


pdf('/home/whou10/data/whou/metpred/evaluate/pbmc/predsc/imp/plot/mefit_hm.pdf', 
    width = 12, height = 100)
pheatmap(mefitsig[row_clustering, ], cluster_rows = F, cluster_cols = F, show_rownames = T, show_colnames = F, scale = 'row')
dev.off()

saveRDS(mefitsig[row_clustering, ], '/home/whou10/data/whou/metpred/evaluate/pbmc/predsc/imp/plot/mefit_hm.rds')
## ====================================
## plot promoter expr and DNAM heatmap
## ====================================
c = corfunc(gefitsig, mefitsig)
summary(c)
posid = which(c > 0.7) # 105 genes
negid = which(c < -0.7) # 167 genes

gefitsig.pos = gefitsig[posid,]
gefitsig.neg = gefitsig[negid,]
mefitsig.pos = mefitsig[posid,]
mefitsig.neg = mefitsig[negid,]


## seperate positive and negative 
pheatmap_result <- pheatmap(gefitsig.pos, cluster_rows = T, cluster_cols = F, show_rownames = F, show_colnames = F, scale = 'row', treeheight_row = 0)
row_clustering <- pheatmap_result$tree_row$order
dev.off()
pdf('/home/whou10/data/whou/metpred/evaluate/pbmc/predsc/imp/plot/gefit_hm_pos.pdf', 
    width = 3, height = 3.8)
pheatmap_result
dev.off()

pdf('/home/whou10/data/whou/metpred/evaluate/pbmc/predsc/imp/plot/mefit_hm_pos.pdf', 
    width = 3, height = 3.8)
pheatmap(mefitsig.pos[row_clustering, ], cluster_rows = F, cluster_cols = F, show_rownames = F, show_colnames = F, scale = 'row')
dev.off()

saveRDS(mefitsig.pos[row_clustering, ], '/home/whou10/data/whou/metpred/evaluate/pbmc/predsc/imp/plot/mefit_hm_pos.rds')

pheatmap_result <- pheatmap(gefitsig.neg, cluster_rows = T, cluster_cols = F, show_rownames = F, show_colnames = F, scale = 'row', treeheight_row = 0)
row_clustering <- pheatmap_result$tree_row$order
dev.off()
pdf('/home/whou10/data/whou/metpred/evaluate/pbmc/predsc/imp/plot/gefit_hm_neg.pdf', 
    width = 3, height = 3.8)
pheatmap_result
dev.off()

pdf('/home/whou10/data/whou/metpred/evaluate/pbmc/predsc/imp/plot/mefit_hm_neg.pdf', 
    width = 3, height = 3.8)
pheatmap(mefitsig.neg[row_clustering, ], cluster_rows = F, cluster_cols = F, show_rownames = F, show_colnames = F, scale = 'row')
dev.off()
saveRDS(mefitsig.neg[row_clustering, ], '/home/whou10/data/whou/metpred/evaluate/pbmc/predsc/imp/plot/mefit_hm_neg.rds')
## ==================
## use dtw distance
## ==================

str(dist)
min(dist)

## select the smallest dtw distance genes
names(dist) = names(pro)[unique(o[,1])]
dist.small <- names(sort(dist)[1:100])
summary(sort(dist)[1:100])
plot(c[dist.small])

dist.small.pos = which(c[dist.small] > 0.7) # 45
dist.small.neg = which(c[dist.small] < -0.7) # 94
str(dist.small.pos)
str(dist.small.neg )
gefit.smalld.pos = gefitsig[dist.small.pos,]
gefit.smalld.neg = gefitsig[dist.small.neg,]
mefit.smalld.pos = mefitsig[dist.small.pos,]
mefit.smalld.neg = mefitsig[dist.small.neg,]

## seperate positive and negative 
pheatmap_result <- pheatmap(gefit.smalld.pos, cluster_rows = T, cluster_cols = F, show_rownames = F, show_colnames = F, scale = 'row', treeheight_row = 0)
row_clustering <- pheatmap_result$tree_row$order
dev.off()
pdf('/home/whou10/data/whou/metpred/evaluate/pbmc/predsc/imp/plot/gefit_hm_smalld_pos.pdf', 
    width = 3, height = 3.8)
pheatmap_result
dev.off()

pdf('/home/whou10/data/whou/metpred/evaluate/pbmc/predsc/imp/plot/mefit_hm_smalld_pos.pdf', 
    width = 3, height = 3.8)
pheatmap(mefit.smalld.pos[row_clustering, ], cluster_rows = F, cluster_cols = F, show_rownames = F, show_colnames = F, scale = 'row')
dev.off()

saveRDS(mefit.smalld.pos[row_clustering, ], '/home/whou10/data/whou/metpred/evaluate/pbmc/predsc/imp/plot/mefit_hm_smalld_pos.rds')

pheatmap_result <- pheatmap(gefit.smalld.neg, cluster_rows = T, cluster_cols = F, show_rownames = F, show_colnames = F, scale = 'row', treeheight_row = 0)
row_clustering <- pheatmap_result$tree_row$order
dev.off()
pdf('/home/whou10/data/whou/metpred/evaluate/pbmc/predsc/imp/plot/gefit_hm_smalld_neg.pdf', 
    width = 3, height = 3.8)
pheatmap_result
dev.off()

pdf('/home/whou10/data/whou/metpred/evaluate/pbmc/predsc/imp/plot/mefit_hm_smalld_neg.pdf', 
    width = 3, height = 3.8)
pheatmap(mefit.smalld.neg[row_clustering, ], cluster_rows = F, cluster_cols = F, show_rownames = F, show_colnames = F, scale = 'row')
dev.off()
saveRDS(mefit.smalld.neg[row_clustering, ], '/home/whou10/data/whou/metpred/evaluate/pbmc/predsc/imp/plot/mefit_hm_smalld_neg.rds')
