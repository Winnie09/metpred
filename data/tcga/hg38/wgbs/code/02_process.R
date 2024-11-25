ddir <- '/home/whou10/data/whou/metpred/data/tcga/hg38/wgbs/raw/download/'
rdir <- '/home/whou10/data/whou/metpred/data/tcga/hg38/wgbs/raw/processed/'
af <- list.files(path=ddir,pattern='.bed.gz')
f = af[1]
me <- lapply(af, function(f){
  d = data.table::fread(paste0(ddir, f), data.table=F)
  v = d[,4]
  names(v) <- paste0(d[,1], ':', d[,2], '_', d[,3])
  v
})
saveRDS(me, paste0(rdir, 'me_cpg_by_sample_list.rds'))  

me2 = do.call(cbind, me)
colnames(me2) = sub('.bed.gz', '',af)
saveRDS(me2, paste0(rdir, 'me_cpg_by_sample.rds'))  

