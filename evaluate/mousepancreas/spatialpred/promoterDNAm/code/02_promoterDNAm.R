rm(list=ls())
ddir = '/home/whou10/data/whou/metpred/evaluate/mousepancreas/spatialpred/imp/res/'
as = list.files(ddir)
id = as.numeric(commandArgs(trailingOnly = T)[[1]])
s = as[id]
af = list.files(paste0(ddir, s))
af = af[!grepl('mean', af)]

#####################################################################
setwd('/home/whou10/data/whou/metpred/')
df = readRDS('/home/whou10/data/whou/metpred/evaluate/mousepancreas/spatialpred/promoterDNAm/res/df_of_cpg_gene.rds')

################################################
rdir = '/home/whou10/data/whou/metpred/evaluate/mousepancreas/spatialpred/promoterDNAm/res/'
d = lapply(af, function(f){
  print(f)
  tmp = readRDS(paste0(ddir, s, '/', f))
  tmp = tmp[rownames(tmp)%in%df[,1], ]
})
d = do.call(rbind, d)
ag = df[,2]
names(ag) = df[,1]
ag = ag[rownames(d)]
rs = rowsum(d, ag)
tab = table(ag)
tab = tab[rownames(rs)]
res = rs/as.vector(tab)
saveRDS(res, paste0(rdir, s, '.rds'))

rm(list=ls())
