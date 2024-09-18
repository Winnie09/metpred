rdir = '/home/whou10/data/whou/metpred/evaluate/pbmc/predsc/imp/res/pred/allcpg/'
af = list.files(rdir, pattern  = 'pred_cpggroup')

source('/home/whou10/scratch16/whou10/resource/myfunc/01_function.R')
ddir <- '/home/whou10/data/whou/metpred/data/pbmc/imp/res/'
rdir2 = '/home/whou10/data/whou/metpred/evaluate/pbmc/predsc/imp/res/pb/'
w.te = readRDS(paste0(ddir, 'wgbs_test.rds'))

for (i in af){
  print(i)
  pred = readRDS(paste0(rdir, i))
  print(str(pred))
  w.te2 = w.te[rownames(pred),]
  str(w.te2)
  pred.ct = sub(':.*', '', colnames(pred))
  pred.ct[grepl('_t', pred.ct)] = 't_cells'
  pred.agg = aggregatefunc2(d = pred, by = pred.ct, fun = 'mean')
  str(pred.agg)
  saveRDS(pred.agg, paste0(rdir2, i))
}


## read in the aggregated pseudobulk of the predicted wgbs
## rbind all CpGs
af = list.files(rdir, pattern  = 'pred_cpggroup')
af = af[order(as.numeric(sub('.rds', '',sub('pred_cpggroup','',af))))]
d <- lapply(af, function(i){
  tmp = readRDS(paste0(rdir2, i))
})
pred.agg = do.call(rbind,d)
saveRDS(pred.agg, paste0(rdir2, 'pred_allcpg.rds'))

