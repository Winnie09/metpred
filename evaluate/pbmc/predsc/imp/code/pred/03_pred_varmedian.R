source('/home/whou10/data/whou/metpred/software/trainpredict.R')
source('/home/whou10/scratch16/whou10/resource/startup.R')
ddir <- '/home/whou10/data/whou/metpred/data/pbmc/imp/res/'
w.tr = readRDS(paste0(ddir, 'wgbs_train.rds'))
# > str(w.tr)
# num [1:30673446, 1:66] 0.744 0.743 0.741 0.741 0.74 ...
# - attr(*, "dimnames")=List of 2
# ..$ : chr [1:30673446] "chr1_10469" "chr1_10471" "chr1_10484" "chr1_10489" ...
# ..$ : chr [1:66] "blueprint_venous_blood-C000S5-CD14-positive_CD16-negative_classical_monocyte" "blueprint_venous_blood-C000S5-mature_neutrophil" "blueprint_venous_blood-C0010K-CD14-positive_CD16-negative_classical_monocyte" "blueprint_venous_blood-C0010K-mature_neutrophil" ...
r.tr = readRDS(paste0(ddir, 'rna_train.rds'))
# > str(r.tr)
# num [1:58434, 1:66] 0 0 2.087 0.566 0.506 ...
# - attr(*, "dimnames")=List of 2
# ..$ : chr [1:58434] "TSPAN6" "TNMD" "DPM1" "SCYL3" ...
# ..$ : chr [1:66] "blueprint_venous_blood-C000S5-CD14-positive_CD16-negative_classical_monocyte" "blueprint_venous_blood-C000S5-mature_neutrophil" "blueprint_venous_blood-C0010K-CD14-positive_CD16-negative_classical_monocyte" "blueprint_venous_blood-C0010K-mature_neutrophil" ...
r.te = readRDS(paste0(ddir, 'rna_test_sc.rds'))
# num [1:8209, 1:37740] 0.0963 0.3414 0.0868 0.0426 0.1903 ...
# - attr(*, "dimnames")=List of 2
# ..$ : chr [1:8209] "NOC2L" "ISG15" "TNFRSF18" "TNFRSF4" ...
# ..$ : chr [1:37740] "b_cells:AAACATACAATGCC-1" "b_cells:AAACATACACGCAT-1" "b_cells:AAACATACGAATAG-1" "b_cells:AAACATACGTGTCA-1" ...
# selcpg = readRDS('/home/whou10/data/whou/metpred/evaluate/pbmc/imp/res/selcpg/top1e4_mono_t_diff_greater0.2_tr_te.rds')
ct = sub(':.*', '', colnames(r.te))
r.te = r.te[, ct %in% c('b_cells', 'cd14_monocytes', 'cd56_nk', 
                        'cytotoxic_t', 'memory_t', 'naive_t',
                        'regulatory_t') ]
ct = sub(':.*', '', colnames(r.te))
id <- unlist(sapply(unique(ct), function(uct){
  which(ct == uct)[1:300]
}, simplify = F))
r.te.sub = r.te[,id]

rm = rowMeans(w.tr)
rvar = apply(w.tr, 1, var)

selcpg = names(rvar > median(rvar, na.rm = T))

mean(selcpg %in% rownames(w.tr))

print(str(selcpg))

pred <- trainpredict(trainexpr=r.tr,
                     testexpr=r.te.sub,
                     trainmeth=w.tr[selcpg,],
                     clunumlist = c(1e3),
                     lambdalist = c(10e-1))
saveRDS(pred, paste0('/home/whou10/data/whou/metpred/evaluate/pbmc/predsc/imp/res/pred/pred_varmedian_', length(selcpg),'_CpG.rds'))