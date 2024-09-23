ddir = '/home/whou10/data/whou/metpred/data/pbmc/imp/res/'
r.te = readRDS(paste0(ddir, 'rna_test_sc.rds'))
str(r.te)
# num [1:8209, 1:37740] 0.0963 0.3414 0.0868 0.0426 0.1903 ...
# - attr(*, "dimnames")=List of 2
# ..$ : chr [1:8209] "NOC2L" "ISG15" "TNFRSF18" "TNFRSF4" ...
# ..$ : chr [1:37740] "b_cells:AAACATACAATGCC-1" "b_cells:AAACATACACGCAT-1" "b_cells:AAACATACGAATAG-1" "b_cells:AAACATACGTGTCA-1" ...
# selcpg = readRDS('/home/whou10/data/whou/metpred/evaluate/pbmc/imp/res/selcpg/top1e4_mono_t_diff_greater0.2_tr_te.rds')
ct = sub(':.*', '', colnames(r.te))
# > table(ct)
# ct
# b_cells cd14_monocytes        cd56_nk 
# 4033            498           7555 
# cytotoxic_t       memory_t        naive_t 
# 7631           6969           4569 
# regulatory_t 
# 6485 

r.te = r.te[, ct %in% c('b_cells', 'cd14_monocytes', 'cd56_nk', 
                        'cytotoxic_t', 'memory_t', 'naive_t',
                        'regulatory_t') ]
# > table(ct)
# ct
# b_cells cd14_monocytes        cd56_nk 
# 4033            498           7555 
# cytotoxic_t       memory_t        naive_t 
# 7631           6969           4569 
# regulatory_t 
# 6485 
ct = sub(':.*', '', colnames(r.te))

id <- unlist(sapply(unique(ct), function(uct){
  which(ct == uct)[1:300]
}, simplify = F))
r.te.sub = r.te[,id]
# > str(r.te.sub)
# num [1:8209, 1:2100] 0.0963 0.3414 0.0868 0.0426 0.1903 ...
# - attr(*, "dimnames")=List of 2
# ..$ : chr [1:8209] "NOC2L" "ISG15" "TNFRSF18" "TNFRSF4" ...
# ..$ : chr [1:2100] "b_cells:AAACATACAATGCC-1" "b_cells:AAACATACACGCAT-1" "b_cells:AAACATACGAATAG-1" "b_cells:AAACATACGTGTCA-1" ...
saveRDS(r.te.sub, paste0(ddir, 'rna_test_sc_sub300.rds'))
