d2 = readRDS('/home/whou10/data/whou/metpred/data/tcga/hg38/2024/combine/EPIC.rds')
write.csv(d2, quote = T, file = '/home/whou10/data/whou/metpred/data/tcga/hg38/2024/combine/EPIC.csv')
str(d2)
write.table(rownames(d2), file = '/home/whou10/data/whou/metpred/data/tcga/hg38/2024/combine/epic_cpgnames.txt', quote = F,
            col.names = F, row.names = F)
write.table(colnames(d2), file = '/home/whou10/data/whou/metpred/data/tcga/hg38/2024/combine/epic_samples.txt', quote = F,
            col.names = F, row.names = F)


d1 = readRDS('/home/whou10/data/whou/metpred/data/tcga/hg38/2024/combine/450k.rds')
write.csv(d1, quote = T, file = '/home/whou10/data/whou/metpred/data/tcga/hg38/2024/combine/450k.csv')
str(d1)
write.table(rownames(d1), file = '/home/whou10/data/whou/metpred/data/tcga/hg38/2024/combine/450k_cpgnames.txt', quote = F,
            col.names = F, row.names = F)
write.table(colnames(d1), file = '/home/whou10/data/whou/metpred/data/tcga/hg38/2024/combine/450k_samples.txt', quote = F,
            col.names = F, row.names = F)



d3 = readRDS('/home/whou10/data/whou/metpred/data/tcga/hg38/2024/combine/wgbs.rds')
write.csv(d3, quote = T, file = '/home/whou10/data/whou/metpred/data/tcga/hg38/2024/combine/wgbs.csv')
str(d3)
write.table(rownames(d3), file = '/home/whou10/data/whou/metpred/data/tcga/hg38/2024/combine/wgbs_cpgnames.txt', quote = F,
            col.names = F, row.names = F)
write.table(colnames(d3), file = '/home/whou10/data/whou/metpred/data/tcga/hg38/2024/combine/wgbs_samples.txt', quote = F,
            col.names = F, row.names = F)



ge = readRDS('/home/whou10/data/whou/metpred/data/tcga/hg38/2024/combine/ge.rds')
write.csv(ge, quote = T, file = '/home/whou10/data/whou/metpred/data/tcga/hg38/2024/combine/ge.csv')
str(ge)
identical(colnames(ge), c(colnames(d1), colnames(d2), colnames(d3)))







