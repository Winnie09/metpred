setwd('/home/whou10/data/whou/metpred/')
rdir <- paste0('evaluate/mousepancreas/trainDuctAcinar_predADM/res/')
rdir.gs <- paste0('evaluate/mousepancreas/trainDuctAcinar_predADM/gs/')
dir.create(rdir.gs, showWarnings = F, recursive = T)

dir.create(rdir, showWarnings = F, recursive = T)
r <- readRDS('data/mousepancreas/spatial/proc/all/cse.rds')
w <- readRDS('data/mousepancreas/wgbs/bs.rds')
int <- intersect(colnames(r),colnames(w))
samp <- sub(':.*','',int)
testid <- int[grep('adm',int)]
trainid <- setdiff(int, testid)
source('software/trainpredict.R')
source('/home/whou10/scratch16/whou10/resource/startup.R')

## group the CpGs from its standard deviation in testing samples smallest to largets slots
## run the prediction for each CpG group
rs <- rowsds(w[,testid])
rs2 <- cut(rs,seq(0,1,0.05))
names(rs2) <- names(rs)

for (rs2.s in unique(rs2)){
  print(rs2.s)
  fn <- paste0(rdir,'/testsd_cpg_', sub('\\]', '', sub('\\(','',rs2.s)), '.rds')
  if (file.exists(fn)){next} 
  pred <- trainpredict(trainexpr=r[,trainid],testexpr=r[,testid],trainmeth=w[names(rs2)[which(rs2 == rs2.s)], trainid])
  saveRDS(pred,file=fn)
  # file = sub('res','gs', fn))
}

