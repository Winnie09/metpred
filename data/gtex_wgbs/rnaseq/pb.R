s <- data.table::fread('GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt',data.table=F)
af <- list.files('/hpc/group/jilab/zj/metpred/data/gtex_wgbs/wgbs/bw')
af <- sub('.mCG.small_smooth.bw','',af)
af2 <- sub('.rds','',list.files('/hpc/group/jilab/zj/resource/eqtl/V7/proc/res'))
af <- intersect(af,af2)
tissue <- s$SMTSD 
tissue <- gsub('_$','',gsub('__','_',gsub(' |\\(|\\)','_',gsub(' - ','_',tissue))))
d <- readRDS('/hpc/group/jilab/zj/metpred/data/gtex_wgbs/rnaseq/log2tpm.rds')
mean(colnames(d)%in%s$SAMPID)

pd <- sapply(af,function(f) {
  print(f)
  sid <- intersect(colnames(d),s[tissue==f,'SAMPID'])
  rowMeans(d[,sid])
})

saveRDS(pd,file='/hpc/group/jilab/zj/metpred/data/gtex_wgbs/rnaseq/pb.rds')
