rm(list=ls())
rdir = '/home/whou10/data/whou/metpred/data/pbmc/imp/res/'
dir.create(rdir, recursive = T)
r = readRDS('/home/whou10/data/whou/metpred/data/v_20220624/combine/rna/hg38_blueprint.rds')
r = r[!is.na(rownames(r)), ]
cn = colnames(r)
cn = cn[!grepl('Leukemia', cn) & !grepl('Myeloma', cn) & !grepl('Lymphoma', cn)]
r = r[, cn]
# pbmc2 = readRDS('/home/whou10/data/whou/rna_imputation/data/processed/pbmc/sorted/norm_genebycell.rds')
pbmc = readRDS('/home/whou10/data/whou/rna_imputation/result/procimpute/pbmc/saver/sorted.rds')
ct = sub(':.*','',colnames(pbmc))
testct <- c('natural_killer_cell',
            'classical_monocyte',
            'B_cell',
            'CD4-positive_alpha-beta_T_cell',
            'CD8-positive_alpha-beta_T_cell')

testid <- trainid <- NULL
for (i in testct){
  tmp = grep(i,cn)
  trainid = c(trainid, tmp[1:ceiling(length(tmp)/2)])
  if (ceiling(length(tmp)/2)+1 != length(tmp)){
    testid = c(testid, tmp[(ceiling(length(tmp)/2)+1) : length(tmp)])
  } else {
    testid = c(testid, tmp[length(tmp)])
  }
}
trainid = setdiff(1:length(cn), testid)
r.tr = r[, trainid]
r.te = r[,testid]
# w = readRDS('/home/whou10/data/whou/metpred/data/v_20220624/combine/wgbs/hg38_blueprint.rds')
# w = readRDS('/home/whou10/data/whou/metpred/data/v_20220624/hg38/blueprint/proc/wgbs_bs.rds')
w = readRDS('/home/whou10/data/whou/metpred/data/v_20220624/hg38/blueprint/proc/wgbs_bs.rds')
w.tr = w[, colnames(r.tr)]
w.te = w[, colnames(r.te)]

## remove cord blood samples in testing gs b/c diff from pbmc 
r.te = r.te[, -c(1,3,4,9)]
w.te = w.te[, -c(1,3,4,9)]
sc.te = pbmc[, ct %in% c('cd56_nk','cd14_monocytes', 'b_cells', 'cytotoxic_t', 'memory_t', 'regulatory_t', 'naive_t')]


saveRDS(r.tr, paste0(rdir, 'rna_train.rds'))
saveRDS(w.tr, paste0(rdir, 'wgbs_train.rds'))
saveRDS(r.te, paste0(rdir, 'rna_test_bulk.rds'))
saveRDS(w.te, paste0(rdir, 'wgbs_test.rds'))
saveRDS(sc.te, paste0(rdir, 'rna_test_sc.rds'))

meta = data.frame(Sample = c(colnames(r.tr), colnames(r.te)), 
                  Usage = c(rep('Train', ncol(r.tr)), rep('Test_goldstandard', ncol(r.te))))
write.csv(meta, '/home/whou10/data/whou/metpred/data/pbmc/imp/res/sampleInfo.csv', row.names = F)

