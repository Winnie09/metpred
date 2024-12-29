## check tcga gene expression batch effect
## across those corresponding to three DNA methylation tech
rdir = '/home/whou10/data/whou/metpred/data/tcga/hg38/2024/quality/res/'
pdir = '/home/whou10/data/whou/metpred/data/tcga/hg38/2024/quality/plot/'

ge = readRDS('/home/whou10/data/whou/metpred/data/tcga/hg38/2024/combine/ge.rds')
wgbssp = readLines('/home/whou10/data/whou/metpred/data/tcga/hg38/2024/combine/wgbs_samples.txt')
epicsp = readLines('/home/whou10/data/whou/metpred/data/tcga/hg38/2024/combine/epic_samples.txt')
arraysp = readLines('/home/whou10/data/whou/metpred/data/tcga/hg38/2024/combine/450k_samples.txt')

ge.wgbs = ge[, wgbssp]
ge.epic = ge[, epicsp]
ge.array = ge[, arraysp]

saveRDS(ge.wgbs, '/home/whou10/data/whou/metpred/data/tcga/hg38/2024/combine/ge_for_wgbs.rds')
saveRDS(ge.epic, '/home/whou10/data/whou/metpred/data/tcga/hg38/2024/combine/ge_for_epic.rds')
saveRDS(ge.array, '/home/whou10/data/whou/metpred/data/tcga/hg38/2024/combine/ge_for_450k.rds')

write.csv(ge.wgbs, quote = F, file = '/home/whou10/data/whou/metpred/data/tcga/hg38/2024/combine/ge_for_wgbs.csv')
write.csv(ge.epic, quote = F, file = '/home/whou10/data/whou/metpred/data/tcga/hg38/2024/combine/ge_for_epic.csv')
write.csv(ge.array, quote = F, file = '/home/whou10/data/whou/metpred/data/tcga/hg38/2024/combine/ge_for_450k.csv')
