library(data.table)
s <- fread('GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt',data.table=F)
saveRDS(s,file='sample.rds')
