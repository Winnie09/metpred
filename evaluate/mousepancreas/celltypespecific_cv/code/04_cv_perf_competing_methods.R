ddir <- '/home/whou10/data/whou/metpred/evaluate/mousepancreas/spatialpred/imp/res/'
meta = read.csv('/home/whou10/data/whou/FeinbergLab/mousePancreas/wgbs/doc/all_adm_sample_metadata.csv')
str(meta)

meta[,3]

w <- readRDS('data/mousepancreas/wgbs/bs.rds')
pred = readRDS('/home/whou10/data/whou/metpred/evaluate/mousepancreas/celltypespecific_cv/res/A46_F_P_K4_D4/testsd_cpg_13.rds')

str(pred) ## num [1:955541, 1:8]
w <- readRDS('data/mousepancreas/wgbs/bs.rds')
str(w)
gs <- w[, colnames(pred)] # 
gs <- gs[rownames(pred), ] # [1:955541, 1:8]
source('/home/whou10/scratch16/whou10/resource/startup.R')
str(gs)



## gs: measured
## predicted: cell-type-specific nnls-based results

## caerulein, normal, D2, 4, 7 [acinar, duct]
## acinar: D4 - D2, D7 - D4
## duct: D4 - D2, D7 - D4

s1 = meta[meta[,2] == 'Caerulein_ADM' & meta[,'Phenotype'] == 'Normal' & meta[,3] == 'Day 2', 1]
s2 = meta[meta[,2] == 'Caerulein_ADM' & meta[,'Phenotype'] == 'Normal' & meta[,3] == 'Day 4', 1]
s3 = meta[meta[,2] == 'Caerulein_ADM' & meta[,'Phenotype'] == 'Normal' & meta[,3] == 'Day 7', 1]
s1a = paste0(s1,':','acinar')
s2a = paste0(s2,':','acinar')
s3a = paste0(s3,':','acinar')


s1, s2, s3
allf = list.files((paste0('/home/whou10/data/whou/metpred/evaluate/mousepancreas/celltypespecific_cv/res/', s1[1])))
allf
pred <- lapply(allf, function(f){
  tmp = readRDS(paste0('/home/whou10/data/whou/metpred/evaluate/mousepancreas/celltypespecific_cv/res/', s1[1], '/',f))
})
str(pred)
pred = do.call(rbind, pred)


allf1.2 = list.files((paste0('/home/whou10/data/whou/metpred/evaluate/mousepancreas/celltypespecific_cv/res/', s1[2])))
allf1.2
pred1.2 <- lapply(allf1.2, function(f){
  tmp = readRDS(paste0('/home/whou10/data/whou/metpred/evaluate/mousepancreas/celltypespecific_cv/res/', s1[2], '/',f))
})
str(pred1.2)
pred1.2 = do.call(rbind, pred1.2)

allf2 = list.files((paste0('/home/whou10/data/whou/metpred/evaluate/mousepancreas/celltypespecific_cv/res/', s2)))
allf
pred2 <- lapply(allf2, function(f){
  tmp = readRDS(paste0('/home/whou10/data/whou/metpred/evaluate/mousepancreas/celltypespecific_cv/res/', s2, '/',f))
})
str(pred2)
pred2 = do.call(rbind, pred2)

allf3 = list.files((paste0('/home/whou10/data/whou/metpred/evaluate/mousepancreas/celltypespecific_cv/res/', s3)))
allf3 ## check job completeness
pred3 <- lapply(allf3, function(f){
  tmp = readRDS(paste0('/home/whou10/data/whou/metpred/evaluate/mousepancreas/celltypespecific_cv/res/', s3, '/',f))
})
pred3 = do.call(rbind, pred3)
str(pred3)
int = Reduce(intersect, list(rownames(pred), rownames(pred1.2), rownames(pred2), rownames(pred3)))
str(int)

pred.all = cbind(pred[int, ,drop=F], pred1.2[int,,drop=F], pred2[int,,drop=F], pred3[int,,drop =F])
str(pred.all)
ct = sub('.*:','',colnames(pred.all))
pred.acinar.diff = pred.all[,ct=='acinar'] %*% matrix(c(-1/2, -1/2, 1,0,0,0, -1,1), ncol = 2)
pred.duct.diff = pred.all[,ct=='duct'] %*% matrix(c(-1/2, -1/2, 1,0,0,0, -1,1), ncol = 2)
str(pred.acinar.diff)

int = intersect(rownames(w), rownames(pred.all))
str(int)
gs.acinar <- w[int, paste0(c(s1,s2,s3),':acinar')] %*% matrix(c(-1/2, -1/2, 1,0,0,0, -1,1), ncol = 2)
gs.duct <- w[int, paste0(c(s1,s2,s3),':duct')] %*% matrix(c(-1/2, -1/2, 1,0,0,0, -1,1), ncol = 2)


pred.acinar.diff = pred.acinar.diff[int,]
pred.duct.diff = pred.duct.diff[int,]

colnames(pred.acinar.diff) <- colnames(gs.acinar) <- paste0('Caerulein:Normal:Acinar:',c('D4-D2','D7-D4'))
colnames(pred.duct.diff) <- colnames(gs.duct) <- paste0('Caerulein:Normal:Duct:',c('D4-D2','D7-D4'))
str(gs.acinar)
a = corfunc[t(pred.acinar.diff), t(gs.acinar)]

i = 1
cor(pred.acinar.diff[,i], gs.acinar[,i])
cor(pred.duct.diff[,i], gs.duct[,i])
i = 2
cor(pred.acinar.diff[,i], gs.acinar[,i])
cor(pred.duct.diff[,i], gs.duct[,i])





## caerulein, pancreatitis, D2, 4, 7 [D2/4, acinar/duct; D7, ADM/duct]
## duct: D4 - D2, D7 - D4


## KLF4, normal, D2,4,7 [acinar, duct]
## acinar: D4 - D2, D7 - D4
## duct: D4 - D2, D7 - D4


## KLF4, pancreatitis, D2, 4, 7 [D2/4, acinar/duct; D7, ADM/duct]
## duct: D4 - D2, D7 - D4

