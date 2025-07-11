library(data.table)
d1 <- fread('/home/whou10/data/whou/metpred/data/tcga/hg38/2024/meta/array/gdc_sample_sheet.2024-12-22.tsv',data.table=F)
d2 <- fread('/home/whou10/data/whou/metpred/data/tcga/hg38/2024/meta/ge/gdc_sample_sheet.2024-12-22.tsv',data.table=F)
d3 <- fread('/home/whou10/data/whou/metpred/data/tcga/hg38/2024/meta/wgbs/tcga_wgbs_meta.csv',data.table=F)

arraycase <- d1[,'Case ID']
gecase <- d2[,'Case ID']
wgbscase <- sub('\t','',setdiff(d3[,2],'not found'))
mean(wgbscase %in% arraycase)

filtercase <- intersect(gecase,arraycase)

d1 <- d1[d1[,'Case ID'] %in% filtercase,]
d2 <- d2[d2[,'Case ID'] %in% filtercase,]

saveRDS(list(array=d1,ge=d2),file='/home/whou10/data/whou/metpred/data/tcga/hg38/2024/meta/combine/filtercase.rds')

