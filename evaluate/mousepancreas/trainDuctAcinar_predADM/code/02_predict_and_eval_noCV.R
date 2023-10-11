setwd('/home/whou10/data/whou/metpred/')
rdir <- paste0('evaluate/mousepancreas/trainDuctAcinar_predADM/res_clu1e3_lambda10e-1/')
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
  fn <- paste0(rdir,'testsd_cpg_', sub('\\]', '', sub('\\(','',rs2.s)), '.rds')
  if (file.exists(fn)){next} 
  pred <- trainpredict(trainexpr=r[,trainid],testexpr=r[,testid],trainmeth=w[names(rs2)[which(rs2 == rs2.s)], trainid], 
                       clunumlist = 1000,
                       lambdalist = 10e-1)
  saveRDS(pred,file=fn)
  # file = sub('res','gs', fn))
}

## evaluate
setwd('/home/whou10/data/whou/metpred/')
w <- readRDS('data/mousepancreas/wgbs/bs.rds')
source('/home/whou10/scratch16/whou10/resource/startup.R')
for (fn in list.files('evaluate/mousepancreas/trainDuctAcinar_predADM/res_clu1e3_lambda10e-1/')){
  print(fn)
  rdir1 = 'evaluate/mousepancreas/trainDuctAcinar_predADM/perf_clu1e3_lambda10e-1/pcc/'
  rdir2 = 'evaluate/mousepancreas/trainDuctAcinar_predADM/perf_clu1e3_lambda10e-1/mse/'
  dir.create(rdir1, recursive = T)
  dir.create(rdir2, recursive = T)
  
  pred <- readRDS(paste0('evaluate/mousepancreas/trainDuctAcinar_predADM/res_clu1e3_lambda10e-1/', fn))
  str(pred) 
  
  gs <- w[, colnames(pred)]  
  gs <- gs[rownames(pred), ] 
  
  str(gs)
  cc <- corfunc(pred, gs)  
  summary(cc)  
  saveRDS(cc, paste0(rdir1, '/', fn))
  
  mse <- sapply(1:nrow(pred), function(i){
    tmpv <- pred[i, ] - gs[i, ]
    mean(tmpv * tmpv)
  })
  names(mse) <- rownames(pred)
  saveRDS(mse, paste0(rdir2, '/', fn))
}


### plot PCC
setwd('/home/whou10/data/whou/metpred/')
source('/home/whou10/scratch16/whou10/resource/startup.R')
dir <- 'evaluate/mousepancreas/trainDuctAcinar_predADM/perf_clu1e3_lambda10e-1/pcc/'
allf <- list.files(dir)
cc <- lapply(allf, function(f){
  print(f)
  tmp <- readRDS(paste0(dir, f))
  chr <- sub(':.*', '', sub('_.*', '', names(tmp)))
  tmp <- tmp[chr %in% paste0('chr', c(seq(1,22), 'X', 'Y'))]
  set.seed(12345)
  tmp <- sample(tmp, min(1e2, length(tmp)), replace = F)
  tmp2 <- data.frame(cor = tmp,  
                     sd = paste0('[',sub('.rds', '', sub('testsd_cpg_', '', f)),')'))
})
cc2 <- do.call(rbind, cc)
acrossRowCor_plot(plotdata = cc2,
                  ylab = 'across-sample PCC',
                  savefile = TRUE,
                  filename = 'evaluate/mousepancreas/trainDuctAcinar_predADM/plot/across_sample_pcc_clu1e3_lambda10e-1.pdf')


### plot MSE
setwd('/home/whou10/data/whou/metpred/')
source('/home/whou10/scratch16/whou10/resource/startup.R')
dir <- 'evaluate/mousepancreas/trainDuctAcinar_predADM/perf_clu1e3_lambda10e-1/mse/'
allf <- list.files(dir)
cc <- lapply(allf, function(f){
  tmp <- readRDS(paste0(dir, f))
  chr <- sub(':.*', '', sub('_.*', '', names(tmp)))
  tmp <- tmp[chr %in% paste0('chr', c(seq(1,22), 'X', 'Y'))]
  set.seed(12345)
  tmp <- sample(tmp, min(1e2, length(tmp)), replace = F)
  tmp2 <- data.frame(cor = tmp,  
                     sd = paste0('[',sub('.rds', '', sub('testsd_cpg_', '', f)),')'))
})
cc2 <- do.call(rbind, cc)
acrossRowCor_plot(plotdata = cc2,
                  ylab = 'across-sample MSE',
                  savefile = TRUE,
                  filename = 'evaluate/mousepancreas/trainDuctAcinar_predADM/plot/across_sample_mse_clu1e3_lambda10e-1.pdf')


