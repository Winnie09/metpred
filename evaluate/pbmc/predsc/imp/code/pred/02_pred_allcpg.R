cpggroup = as.numeric(commandArgs(trailingOnly = T)[1])
source('/home/whou10/data/whou/metpred/software/trainpredict.R')
source('/home/whou10/scratch16/whou10/resource/startup.R')
ddir <- '/home/whou10/data/whou/metpred/data/pbmc/imp/res/'
w.tr = readRDS(paste0(ddir, 'wgbs_train.rds'))
str(w.tr)
# num [1:30673446, 1:66] 0.744 0.743 0.741 0.741 0.74 ...
# - attr(*, "dimnames")=List of 2
# ..$ : chr [1:30673446] "chr1_10469" "chr1_10471" "chr1_10484" "chr1_10489" ...
# ..$ : chr [1:66] "blueprint_venous_blood-C000S5-CD14-positive_CD16-negative_classical_monocyte" "blueprint_venous_blood-C000S5-mature_neutrophil" "blueprint_venous_blood-C0010K-CD14-positive_CD16-negative_classical_monocyte" "blueprint_venous_blood-C0010K-mature_neutrophil" ...
rs2 = readRDS(paste0(ddir, 'wgbs_train_nonNA_CpG_rowsds_cut.rds'))
w.tr = w.tr[names(rs2)[rs2 == levels(rs2)[cpggroup]], , drop = F]
if (nrow(w.tr) < 1){
  stop('No CpG in this rowsds group!')
}

r.tr = readRDS(paste0(ddir, 'rna_train.rds'))
str(r.tr)
# num [1:58434, 1:66] 0 0 2.087 0.566 0.506 ...
# - attr(*, "dimnames")=List of 2
# ..$ : chr [1:58434] "TSPAN6" "TNMD" "DPM1" "SCYL3" ...
# ..$ : chr [1:66] "blueprint_venous_blood-C000S5-CD14-positive_CD16-negative_classical_monocyte" "blueprint_venous_blood-C000S5-mature_neutrophil" "blueprint_venous_blood-C0010K-CD14-positive_CD16-negative_classical_monocyte" "blueprint_venous_blood-C0010K-mature_neutrophil" ...
r.te.sub = readRDS(paste0(ddir, 'rna_test_sc_sub300.rds'))
 
pred <- trainpredict(trainexpr=r.tr,
                     testexpr=r.te.sub,
                     trainmeth=w.tr,
                     clunumlist = c(1e3),
                     lambdalist = c(10e-1))
saveRDS(pred, paste0('/home/whou10/data/whou/metpred/evaluate/pbmc/predsc/imp/res/pred/allcpg/pred_cpggroup', cpggroup, '.rds'))
