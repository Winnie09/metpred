# 1. on training me, select cpg with no NA. 
# 2. on no NA cpg, calculate the mandiff between HSC and MEP. save meandiff, further select cpg with largest meandiff 
w = readRDS('/home/whou10/data/whou/metpred/data/hsc/proc/wgbs.rds')
ct = sapply(colnames(w), function(i) strsplit(i, '_')[[1]][2])
w2 = w[, ct %in% c('HSC', 'MPP', 'CMP', 'MEP')]
w.hsc = w[, ct %in% c('HSC')]
w.mep = w[, ct %in% c('MEP')]
diff = rowMeans(w.mep, na.rm = T) - rowMeans(w.hsc, na.rm = T)
saveRDS(diff, '/home/whou10/data/whou/metpred/evaluate/hsc/erythroid/selcpg/mep_hsc/meandiff.rds')

me = readRDS(file='/home/whou10/data/whou/metpred/data/hsc/proc/wgbs_matchedtrain.rds')
me2 = me[complete.cases(me),]
int = intersect(rownames(w), rownames(me2))
diff.abs = abs(diff[int])
selcpg = names(sort(diff.abs, decreasing = T))[1:1e4]
saveRDS(selcpg, '/home/whou10/data/whou/metpred/evaluate/hsc/erythroid/selcpg/mep_hsc/meandiff_top1e4cpg_noNA_in_training.rds')

