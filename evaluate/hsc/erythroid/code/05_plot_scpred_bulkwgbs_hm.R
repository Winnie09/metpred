pred2 = readRDS(file='/home/whou10/data/whou/metpred/evaluate/hsc/erythroid/res/pred_rd1e5.rds')
library(pheatmap)
png('/home/whou10/data/whou/metpred/evaluate/hsc/erythroid/plot/pred_hm.png')
pheatmap(pred2, cluster_rows = F, cluster_cols = F, show_colnames = F, show_rownames = F, labels_col = 'Single cells ordered by pseudotime', labels_row = 'top 1e4 CpGs largest mean difference')
dev.off()

w3 = readRDS('/home/whou10/data/whou/metpred/evaluate/hsc/erythroid/bulk/bulk_wgbs.rds')
library(pheatmap)
png('/home/whou10/data/whou/metpred/evaluate/hsc/erythroid/plot/bulk_wgbs_hm.png')
pheatmap(w3, cluster_rows = F, cluster_cols = F, show_rownames = F, labels_row = 'top 1e4 CpGs largest mean difference')
dev.off()

me = readRDS(file='/home/whou10/data/whou/metpred/data/hsc/proc/wgbs_matchedtrain.rds')
library(pheatmap)
png('/home/whou10/data/whou/metpred/evaluate/hsc/erythroid/plot/wgbs_training_top1e4cpg.png')
pheatmap(me[rownames(pred2), ], cluster_rows = F, cluster_cols = F, show_rownames = F, labels_row = 'top 1e4 CpGs largest mean difference')
dev.off()

