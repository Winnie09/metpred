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

## ======== <<<<<<<<<< filter out chrX and chrY
d1.chr = sub('_.*','', rownames(d1))
d1 = d1[d1.chr %in% paste0('chr', 1:22), ]
print('450k dimension')
print(dim(d1))

d2.chr = sub('_.*','', rownames(d2))
d2 = d2[d2.chr %in% paste0('chr', 1:22), ]
print('EPIC dimension')
print(dim(d2))

d3.chr = sub('_.*','', rownames(d3))
d3 = d3[d3.chr %in% paste0('chr', 1:22), ]
print('WGBS dimension')
print(dim(d3))

## ==================== >>>>>>>>>>>

mean(colnames(ge) %in% c(colnames(d1),colnames(d2),colnames(d3)))
mean(colnames(d1) %in% colnames(ge))
mean(colnames(d2) %in% colnames(ge))
mean(colnames(d3) %in% colnames(ge))

mean(rownames(d1) %in% rownames(d2))
mean(rownames(d1) %in% rownames(d3))
mean(rownames(d2) %in% rownames(d3))


na.mean = colMeans(is.na(d2))
summary(na.mean)
unname(sort(na.mean))
id = which(na.mean <= 0.3)
str(id)
summary(na.mean[id])
d2 = d2[, id ]
str(d2.2)

saveRDS(d1,'/home/whou10/data/whou/metpred/data/tcga/hg38/2024/combine/450k.rds')
saveRDS(d2,'/home/whou10/data/whou/metpred/data/tcga/hg38/2024/combine/EPIC.rds')
saveRDS(d3,'/home/whou10/data/whou/metpred/data/tcga/hg38/2024/combine/wgbs.rds')
# saveRDS(ge,'/home/whou10/data/whou/metpred/data/tcga/hg38/2024/combine/ge.rds')

# 20241224 
ge2 = ge[, colnames(ge) %in% c(colnames(d1), colnames(d2), colnames(d3))]
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




write.table(rownames(d1), file = '/home/whou10/data/whou/metpred/data/tcga/hg38/2024/combine/450k_cpgnames.txt', quote = F,
            col.names = F, row.names = F)
write.table(rownames(d2), file = '/home/whou10/data/whou/metpred/data/tcga/hg38/2024/combine/epic_cpgnames.txt', quote = F,
            col.names = F, row.names = F)
write.table(rownames(d3), file = '/home/whou10/data/whou/metpred/data/tcga/hg38/2024/combine/wgbs_cpgnames.txt', quote = F,
            col.names = F, row.names = F)

write.table(colnames(d1), file = '/home/whou10/data/whou/metpred/data/tcga/hg38/2024/combine/450k_samples.txt', quote = F,
            col.names = F, row.names = F)
write.table(colnames(d2), file = '/home/whou10/data/whou/metpred/data/tcga/hg38/2024/combine/epic_samples.txt', quote = F,
            col.names = F, row.names = F)
write.table(colnames(d3), file = '/home/whou10/data/whou/metpred/data/tcga/hg38/2024/combine/wgbs_samples.txt', quote = F,
            col.names = F, row.names = F)
