d = readRDS('/home/whou10/data/whou/metpred/data/tcga/hg38/wgbs/combine/me_cpg_by_sample.rds')
ct1 = paste0('TCGA-','',sub('_.*','',sub('TCGA_', '', colnames(d))))
df = data.frame(Sample = colnames(d), Cancer = ct1)
write.table(df, '/home/whou10/data/whou/metpred/data/tcga/hg38/wgbs/combine/project.csv', sep = ',', quote = F, row.names = F)

pro = readRDS('/home/whou10/data/whou/metpred/data/tcga/hg38/450k_2021/combine/project.rds')
df2 = data.frame(Sample = names(pro), Cancer = unname(pro))
write.table(df2, '/home/whou10/data/whou/metpred/data/tcga/hg38/450k_2021/combine/Methylformer_data/project.csv', sep = ',', quote = F, row.names = F)


