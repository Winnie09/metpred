## =============
## Prepare TCGA data
## =======================
setwd('/home/whou10/data/whou/metpred/data/tcga/hg38/450k_2021/combine')
loc = readRDS('/home/whou10/data/whou/metpred/data/tcga/hg38/450k_2021/loc/hg38.rds')
me = readRDS('me.rds')
mi = readRDS('meimpute.rds')

# some CpGs have the same location, remove duplicated locations
# cg16566605 - chr1:13410375-13410376
# cg20017108 - chr1:13410375-13410376
# 
# 
# > which(rn == 'chr1_13410375')
# cg16566605 cg20017108 
#     249207     294915 

rn = loc[rownames(me)]      
c = sub('_.*','',rn)
table(c)



me2 = me[c %in% paste0('chr', 1:22), ]
rn = loc[rownames(me2)]      

write.csv(data.frame(CpG_name = rownames(me2), CpG_location = rn, stringsAsFactors = FALSE), 'Methylformer_data/CpG_name_location.csv', row.names = FALSE)

rownames(me2) = rn
write.csv(me2, 'Methylformer_data/me_rownamesloc.csv')
saveRDS(me2, 'Methylformer_data/me_rownamesloc.rds')
saveRDS(rownames(me2), 'Methylformer_data/me_cpgnames.rds')


rn = loc[rownames(mi)]
c = sub('_.*','',rn)
mi2 = mi[c %in% paste0('chr', 1:22), ]
rn2 = loc[rownames(mi2)]
rownames(mi2) = rn2
write.csv(mi2, 'Methylformer_data/meimpute_rownamesloc.csv')
saveRDS(mi2, 'Methylformer_data/meimpute_rownamesloc.rds')


ge = readRDS('ge.rds')
write.csv(ge, 'ge.csv')


