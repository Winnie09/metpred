library(limma)
expr <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/metpred/nullsimu/data/data/4To4_allcell/pred_ery.rds')
expr <- expr[1:10000,]
source('/home-4/zji4@jhu.edu/scratch/raisin/software/raisin/raisin.R')
source('/home-4/whou10@jhu.edu/scratch/Wenpin/metpred/software/findDMR.R')
time <- NULL
func <- func[!names(func)=='limmacell_saver']
rdir <- '/home-4/whou10@jhu.edu/work-zfs/whou10/metpred/scalability/result/'

sn = 5 ## 5, 10 
cn <- as.numeric(commandArgs(trailingOnly = T)[[1]])   # 10,100,1000,10000
print(sn)
print(cn)

sample <- paste0(rep(1:(sn*2),each=cn),'_',rep(1:2,each=sn*cn))
set.seed(12345)
expr <- expr[,sample(1:ncol(expr),length(sample),replace = T)]
colnames(expr) <- paste0('cell',1:ncol(expr))
samplename <- unique(sample)
group <- as.numeric(sub('.*_','',samplename))
ind <- NULL

myfun <- func[['t_saver']]
res <- myfun()
saveRDS(res, paste0(sn, '_', cn,'_res.rds'))

