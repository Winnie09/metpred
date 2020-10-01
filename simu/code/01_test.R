source('/home-4/whou10@jhu.edu/scratch/Wenpin/metpred/software/findDMR.R')
ddir <-  '/home-4/whou10@jhu.edu/scratch/Wenpin/metpred/simu/data/data/'
func <- func[grepl('saver',names(func))]
func <- func[!names(func) %in% 'limmacell_saver']
for (type in rev(list.files(ddir))){
  rdir <- paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/metpred/simu/res/',type,'/')
  if (!file.exists(rdir)){
    print(type)
    expr <- readRDS(paste0('/scratch/users/whou10@jhu.edu/Wenpin/metpred/simu/data/data/', type, '/me.rds'))
    sample <- sub(':.*','',colnames(expr))
    samplename <- unique(sample)
    sum(as.numeric(sub('BM','',samplename)) %in% c(1,2,5,6) + 0 == 1) < length(samplename)
    group <- as.numeric(sub('BM','',samplename)) %in% c(1,2,5,6) + 1
    ind <- NULL
    for (metid in 1:length(func)) {
      fun <- func[[metid]]
      res <- fun()
      dir.create(rdir)
      saveRDS(res,file=paste0(rdir, names(func)[metid],'.rds'))  
    }
    design <- data.frame(sample=samplename,feature=group)
    source('/home-4/whou10@jhu.edu/scratch/Wenpin/raisin/software/raisin/raisin.R')
    res <- RAISINtest(RAISINfit(log(expr),sample,design,testtype='unpaired',filtergene=F))
    saveRDS(res,file=paste0(rdir, 'ourmethod.rds'))
  
  }
}
  
## find number of DMR (all are FP)
rdir <- '/home-4/whou10@jhu.edu/scratch/Wenpin/metpred/simu/res/'
atype <- list.files(rdir)
atype <- atype[atype != 'perf']
num <- sapply(atype, function(type){
    af <- list.files(paste0(rdir, type))
    len.fp <- NULL
    for (f in af){
      res <- readRDS(file=paste0(rdir, type, '/', f))
      len.fp[sub('.rds','',f)] <- nrow(res[res$FDR < 0.05, ])
    }
    return(len.fp)
})
saveRDS(num, '/home-4/whou10@jhu.edu/scratch/Wenpin/metpred/simu/perf/num.fp.rds')



