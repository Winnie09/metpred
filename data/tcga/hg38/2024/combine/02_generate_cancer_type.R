ge = readRDS('/home/whou10/data/whou/metpred/data/tcga/hg38/2024/combine/ge.rds')

pro = read.csv('/home/whou10/data/whou/metpred/data/tcga/hg38/450k_2021/combine/Methylformer_data/project.csv')
prosp = pro$Cancer
names(prosp) = sapply(pro$Sample, function(i) paste0(strsplit(i, '-')[[1]][1:3], collapse = '-') )
prosp = prosp[names(prosp) %in% colnames(ge)]

## previous wgbs sample-cancer table: does not find matched samples
pro2 = read.csv('/home/whou10/data/whou/metpred/data/tcga/hg38/2024/wgbs/combine/project.csv')
prosp2 = pro2$Cancer
names(prosp2) = sapply(pro2$Sample, function(i) paste0(strsplit(i, '_')[[1]][1:3], collapse = '-') )
prosp2 = prosp2[names(prosp2) %in% colnames(ge)]

## include a "cancer" meta
cancer = rep('NA', ncol(ge))
names(cancer) = colnames(ge)
cancer[names(prosp)] = prosp

str(cancer)
tmp = data.frame(Sample = colnames(ge), Cancer = cancer)
str(tmp)
write.csv(tmp, 
          file= '/home/whou10/data/whou/metpred/data/tcga/hg38/2024/combine/project.csv', 
          row.names = F,
          quote = F)
saveRDS(tmp, file= '/home/whou10/data/whou/metpred/data/tcga/hg38/2024/combine/project.rds')

