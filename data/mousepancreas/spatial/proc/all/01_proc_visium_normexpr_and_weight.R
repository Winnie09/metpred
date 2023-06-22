library(Seurat)
library(Matrix)
setwd('/home/whou10/data/whou/FeinbergLab/mousePancreas/visium/deconvolute/RCTD/')
load('adm_spatial_RCTD_manual_and_pseudospot_normalized_weights_matrix.rda')
load('seurat_adm_spatial_for_ji_lab.rda')
meta <- adm.combined@meta.data
meta <- meta[rownames(norm_weights),]

## wgbs and visium samples names not the same, here we match them.
samp <- as.character(meta$sample)
samp[samp=='D2_A29_P_K4'] <- 'A29_M_P_K4_D2'
samp[samp=='D2_A34_P_K4'] <- 'A34_M_P_K4_D2'
samp[samp=='D2_A53_LM'] <- 'A53_F_LM_D2'
samp[samp=='D2_A60_LM'] <- 'A60_M_WT_LM_D2'
samp[samp=='D4_A31_P_K4'] <- 'A31_M_P_K4_D4'
samp[samp=='D4_A46_P_K4'] <- 'A46_F_P_K4_D4'
samp[samp=='D4_A73_LM'] <- 'A73_F_LM_D4'
samp[samp=='D7_A47_P_K4'] <- 'A47_F_P_K4_D7'
samp[samp=='D7_A66_P_K4'] <- 'A66_M_P_K4_D7'
samp[samp=='D7_A72_LM'] <- 'A72_F_LM_D7'

## save the normalized expression,counts and the deconvoluted celltype weights
m <- adm.combined@assays$Spatial@data
m <- m[,rownames(norm_weights)]
cnt <- adm.combined@assays$Spatial@counts
cnt <- cnt[,rownames(norm_weights)]


colnames(cnt) <- colnames(m) <- rownames(norm_weights) <- paste0(samp,':',rownames(norm_weights))
saveRDS(m,file='/home/whou10/data/whou/metpred/data/mousepancreas/spatial/proc/all/normexpr.rds')
saveRDS(norm_weights,file='/home/whou10/data/whou/metpred/data/mousepancreas/spatial/proc/all/weight.rds')
saveRDS(cnt,file='/home/whou10/data/whou/metpred/data/mousepancreas/spatial/proc/all/count.rds')


