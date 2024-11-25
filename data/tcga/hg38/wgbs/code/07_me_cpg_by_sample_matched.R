af = list.files('/home/whou10/data/whou/metpred/data/tcga/hg38/wgbs/raw/download/rna')
af = af[!grepl("\\.sh", af)]
rdir = '/home/whou10/data/whou/metpred/data/tcga/hg38/wgbs/combine/'
ddir = '/home/whou10/data/whou/metpred/data/tcga/hg38/wgbs/processed/'
me = readRDS(paste0(ddir, 'me_cpg_by_sample.rds'))
cn = sub('.bed.gz', '', colnames(me))
colnames(me) = cn
saveRDS(me, paste0(rdir, 'me_cpg_by_sample.rds'))

me2 = me[, colnames(me) %in% af]

saveRDS(me2, paste0(rdir, 'me_cpg_by_sample_matched.rds'))
saveRDS(colnames(me2), paste0(rdir, 'matched_samples.rds'))
write.table(me2, file = paste0(rdir, 'me_cpg_by_sample_matched.txt'), sep = '\t')

