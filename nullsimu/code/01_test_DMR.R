for (type in list.files('/home-4/whou10@jhu.edu/scratch/Wenpin/metpred/nullsimu/data/data/')){
  source('/home-4/whou10@jhu.edu/scratch/Wenpin/metpred/software/findDMR.R')
  expr <- readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/metpred/nullsimu/data/data/', type, '/pred_ery.rds'))
  sample <- sub(':.*','',colnames(expr))
  samplename <- unique(sample)
  if (length(samplename) == 8) group <- as.numeric(sub('BM','',samplename)) %in% c(1,2,5,6) + 1
  if (length(samplename) == 4) group <- as.numeric(sub('BM','',samplename)) %in% c(5,6) + 1
  ind <- NULL
  func <- func[grepl('saver',names(func))]
  
  for (metid in 1:length(func)) {
    print(metid)
    fun <- func[[metid]]
    res <- fun()
    dir.create(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/metpred/nullsimu/res/',type))
    saveRDS(res,file=paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/metpred/nullsimu/res/',type, '/', names(func)[metid],'.rds'))
  }

  design <- data.frame(sample=samplename,feature=group)
  source('/home-4/whou10@jhu.edu/scratch/Wenpin/raisin/software/raisin/raisin.R')
  res <- RAISINtest(RAISINfit(log(expr),sample,design,testtype='unpaired',filtergene=F))
  saveRDS(res,file=paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/metpred/nullsimu/res/', type, '/', 'ourmethod.rds'))
}
  


## find number of DMR (all are FP)
atype <- list.files('/home-4/whou10@jhu.edu/scratch/Wenpin/metpred/nullsimu/res/')
atype <- atype[atype != 'perf']
num <- sapply(atype, function(type){
    af <- list.files(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/metpred/nullsimu/res/', type))
    len.fp <- NULL
    for (f in af){
      res <- readRDS(file=paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/metpred/nullsimu/res/', type, '/', f))
      len.fp[sub('.rds','',f)] <- nrow(res[res$FDR < 0.05, ])
    }
    return(len.fp)
})
saveRDS(num, '/home-4/whou10@jhu.edu/scratch/Wenpin/metpred/nullsimu/res/perf/num.fp.rds')


