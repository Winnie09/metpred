source('/home-4/whou10@jhu.edu/work-zfs/whou10/metpred/software/findDMR.R')
source('/home-4/whou10@jhu.edu/work-zfs/whou10/trajectory_variability/function/01_function.R')
ddir <-  '/home-4/whou10@jhu.edu/work-zfs/whou10/metpred/simu/data/data/'
func <- func[grepl('saver',names(func))]
func <- func[!names(func) %in% 'limmacell_saver']
allf <- list.files(ddir)
getf <- list.files('/home-4/whou10@jhu.edu/work-zfs/whou10/metpred/simu/res/')

if (length(getf) > 0){
  v <- sapply(getf, function(f){
    if (length(list.files(paste0('/home-4/whou10@jhu.edu/work-zfs/whou10/metpred/simu/res/', f))) == 7) return(TRUE)
    else return(FALSE)
  })
  getf <- getf[v]
  runf <- setdiff(allf, getf)
} else {
  runf <- allf
}

runf <- runf[1:20]  

for (type in runf){
  rdir <- paste0('/home-4/whou10@jhu.edu/work-zfs/whou10/metpred/simu/res/',type,'/')
  dir.create(rdir, showWarnings = FALSE)
  if (!file.exists(rdir) | length(list.files(rdir)) != 7){
    print(type)
    expr <- readRDS(paste0('/home-4/whou10@jhu.edu/work-zfs/whou10/metpred/simu/data/data/', type, '/me.rds'))
    sample <- sub(':.*','',colnames(expr))
    samplename <- unique(sample)
    sum(as.numeric(sub('BM','',samplename)) %in% c(1,2,5,6) + 0 == 1) < length(samplename)
    group <- as.numeric(sub('BM','',samplename)) %in% c(1,2,5,6) + 1
    ind <- NULL
    for (metid in 1:length(func)) {
      fun <- func[[metid]]
      res <- fun()
      saveRDS(res,file=paste0(rdir, names(func)[metid],'.rds'))  
    }
    design <- data.frame(sample=samplename,feature=group)
    source('/home-4/whou10@jhu.edu/scratch/Wenpin/raisin/software/raisin/raisin.R')
    res <- RAISINtest(RAISINfit(log(expr),sample,design,testtype='unpaired',filtergene=F))
    saveRDS(res,file=paste0(rdir, 'ourmethod.rds'))
  }
}

