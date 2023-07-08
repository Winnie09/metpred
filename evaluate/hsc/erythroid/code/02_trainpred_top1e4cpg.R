ddir <- '/home/whou10/scratch16/whou10/trajectory_variability/hca_bone_marrow_data_analysis/real/build_from_tree_variability/result/erythroid/'
expr = readRDS(paste0(ddir, 'input_expr.rds'))
rownames(expr) = sub(':.*','',rownames(expr))

ge =  readRDS(file='/home/whou10/data/whou/metpred/data/hsc/proc/rna_matchedtrain.rds')
me = readRDS(file='/home/whou10/data/whou/metpred/data/hsc/proc/wgbs_matchedtrain.rds')

## methylation not BSmooth yet. Remove NA. 
me2 = me[complete.cases(me),]

selcpg = readRDS('/home/whou10/data/whou/metpred/evaluate/hsc/erythroid/selcpg/mep_hsc/meandiff_top1e4cpg_noNA_in_training.rds')


# prediction
source('/home/whou10/data/whou/metpred/software/trainpredict.R')
pred <- trainpredict(trainexpr=ge,testexpr=expr, trainmeth=me2[1:1e5,], clunumlist = 500, lambdalist = 1e-1)   
saveRDS(pred,file='/home/whou10/data/whou/metpred/evaluate/hsc/erythroid/res/pred_rd1e5.rds')



