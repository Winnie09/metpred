source('/home/whou10/scratch16/whou10/resource/myfunc/01_function.R')
ddir = '/home/whou10/data/whou/metpred/evaluate/pbmc/predsc/imp/res/pred/allcpg/'
rdir = '/home/whou10/data/whou/metpred/evaluate/pbmc/predsc/imp/perf/ramp/'
af = list.files(ddir, pattern = 'pred_cpggroup')
f = af[10]
res = lapply(af, function(f){
  print(f)
  d = readRDS(paste0(ddir, f))
  d = d[!is.na(rowSums(d)),]
  rs = rowsds(d)
  rs = rs[-which(is.na(rs))]
  dtmp = d[names(rs)[rs >=0.3], ]
})
res2 = do.call(rbind, res)
saveRDS(res2, paste0(rdir, 'predicted_cpg_sd0.3.rds'))  
