d1 <- readRDS('/home/whou10/data/whou/metpred/data/tcga/hg38/2024/array/mat/case_418858.rds')
p1 <- readRDS('/home/whou10/data/whou/metpred/data/tcga/hg38/2024/array/pos/450k.rds')
d1 <- d1[rownames(d1) %in% names(p1),]
rownames(d1) <- unname(p1[rownames(d1)])


d2 <- readRDS('/home/whou10/data/whou/metpred/data/tcga/hg38/2024/array/mat/case_758023.rds')
p2 <- readRDS('/home/whou10/data/whou/metpred/data/tcga/hg38/2024/array/pos/EPIC.rds')
d2 <- d2[rownames(d2) %in% names(p2),]
rownames(d2) <- unname(p2[rownames(d2)])

d3 <- readRDS('/home/whou10/data/whou/metpred/data/tcga/hg38/2024/wgbs/combine/me_cpg_by_sample.rds')
m <- data.table::fread('/home/whou10/data/whou/metpred/data/tcga/hg38/2024/wgbs/meta/tcga_wgbs_meta.csv',data.table=F)
m <- m[m[,2]!='not found',]
m[,2] <- sub('\t','',m[,2])
rownames(m) <- m[,1]
colnames(d3) <- m[colnames(d3),2]

ge <- readRDS('/home/whou10/data/whou/metpred/data/tcga/hg38/2024/ge/mat/log2tpm_case.rds')

d2 <- d2[,colnames(d2) %in% colnames(ge)]

mean(colnames(ge) %in% c(colnames(d1),colnames(d2),colnames(ge)))
mean(colnames(d1) %in% colnames(ge))
mean(colnames(d2) %in% colnames(ge))
mean(colnames(d3) %in% colnames(ge))

mean(rownames(d1) %in% rownames(d2))
mean(rownames(d1) %in% rownames(d3))
mean(rownames(d2) %in% rownames(d3))

saveRDS(d1,'/home/whou10/data/whou/metpred/data/tcga/hg38/2024/combine/450k.rds')
saveRDS(d2,'/home/whou10/data/whou/metpred/data/tcga/hg38/2024/combine/EPIC.rds')
saveRDS(d3,'/home/whou10/data/whou/metpred/data/tcga/hg38/2024/combine/wgbs.rds')
# saveRDS(ge,'/home/whou10/data/whou/metpred/data/tcga/hg38/2024/combine/ge.rds')

# 20241224 
ge2 = ge[, c(colnames(d1), colnames(d2), colnames(d3))]
rn = rownames(ge2)
gn = sub(':.*','', rn)
gid = sub('.*:','',rn)
table(substr(gid, 1,3))
rownames(ge2) = paste0(gn, ';', gid)
saveRDS(ge2, '/home/whou10/data/whou/metpred/data/tcga/hg38/2024/combine/ge.rds')
write.csv(ge2, quote = F, file = '/home/whou10/data/whou/metpred/data/tcga/hg38/2024/combine/ge.csv')
write.csv(d1, quote = F, file = '/home/whou10/data/whou/metpred/data/tcga/hg38/2024/combine/450k.csv')
write.csv(d2, quote = F, file = '/home/whou10/data/whou/metpred/data/tcga/hg38/2024/combine/EPIC.csv')
write.csv(d3, quote = F, file = '/home/whou10/data/whou/metpred/data/tcga/hg38/2024/combine/wgbs.csv')


