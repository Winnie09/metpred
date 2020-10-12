library(limma)
expr <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/metpred/nullsimu/data/data/4To4_allcell/pred_ery.rds')
expr <- expr[1:10000,]
source('/home-4/zji4@jhu.edu/scratch/raisin/software/raisin/raisin.R')
source('/home-4/whou10@jhu.edu/scratch/Wenpin/metpred/software/findDMR.R')
time <- NULL
func <- func[!names(func)=='limmacell_saver']
rdir <- '/home-4/whou10@jhu.edu/work-zfs/whou10/metpred/scalability/result/'
  
sn <- as.numeric(commandArgs(trailingOnly = T)[[1]])   # 3,4,5,10
print(sn)

for (cn in c(10,100,1000,10000)) {
  print(cn)
  sample <- paste0(rep(1:(sn*2),each=cn),'_',rep(1:2,each=sn*cn))
  set.seed(12345)
  expr <- expr[,sample(1:ncol(expr),length(sample),replace = T)]
  colnames(expr) <- paste0('cell',1:ncol(expr))
  samplename <- unique(sample)
  group <- as.numeric(sub('.*_','',samplename))
  ind <- NULL

  for (metid in names(func)) {
    print(metid)
    myfun <- func[[metid]]
    time_sec <- as.numeric(system.time({myfun()})[3])
    time <- rbind(time,data.frame(met = metid, cn = cn, sn = sn, time_sec = time_sec))
  }
  time_sec <- as.numeric(system.time({RAISINtest(RAISINfit(expr,sample,testtype='unpaired',design=data.frame(sample=samplename,feature=group),ncores=10))})[3])
  time <- rbind(time,data.frame(met = 'ourmethod', cn = cn, sn = sn, time_sec))
}
colnames(time)[4] <- 'time'
saveRDS(time, paste0(rdir, 'time_sn', sn, '.rds'))

