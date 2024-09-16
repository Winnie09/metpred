source('/home/whou10/scratch16/whou10/resource/startup.R')
ddir <- '/home/whou10/data/whou/metpred/data/pbmc/imp/res/'
w.tr = readRDS(paste0(ddir, 'wgbs_train.rds'))
# num [1:30673446, 1:66] 0.744 0.743 0.741 0.741 0.74 ...
# - attr(*, "dimnames")=List of 2
# ..$ : chr [1:30673446] "chr1_10469" "chr1_10471" "chr1_10484" "chr1_10489" ...
# ..$ : chr [1:66] "blueprint_venous_blood-C000S5-CD14-positive_CD16-negative_classical_monocyte" "blueprint_venous_blood-C000S5-mature_neutrophil" "blueprint_venous_blood-C0010K-CD14-positive_CD16-negative_classical_monocyte" "blueprint_venous_blood-C0010K-mature_neutrophil" ...

w.tr = w.tr[complete.cases(w.tr), , drop = F]
# str(w.tr)
# num [1:30672982, 1:66] 0.744 0.743 0.741 0.741 0.74 ...
# - attr(*, "dimnames")=List of 2
# ..$ : chr [1:30672982] "chr1_10469" "chr1_10471" "chr1_10484" "chr1_10489" ...
# ..$ : chr [1:66] "blueprint_venous_blood-C000S5-CD14-positive_CD16-negative_classical_monocyte" "blueprint_venous_blood-C000S5-mature_neutrophil" "blueprint_venous_blood-C0010K-CD14-positive_CD16-negative_classical_monocyte" "blueprint_venous_blood-C0010K-mature_neutrophil" ...
rs = rowsds(w.tr)
rs2 <- cut(rs, seq(0, 1, 0.05))
names(rs2) <- names(rs)
saveRDS(rs2, paste0(ddir, 'wgbs_train_nonNA_CpG_rowsds_cut.rds'))



