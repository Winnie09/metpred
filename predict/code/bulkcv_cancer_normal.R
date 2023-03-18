tb = read.csv('/home-4/whou10@jhu.edu/scratch/Wenpin/metpred/doc/ct_tissue_info_for_matched_wgbs_and_RNAseq.csv', header = TRUE, as.is = TRUE)
source('/home-4/whou10@jhu.edu/scratch/Wenpin/metpred/software/trainmodel.R')
source('/home-4/whou10@jhu.edu/scratch/Wenpin/metpred/software/predict.R')
d.bak <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/metpred/final/procrna.rds')
m.bak <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/metpred/final/nonawgbs_hg38.rds')
rdir <- '/home-4/whou10@jhu.edu/scratch/Wenpin/metpred/predict/res/'

identical(colnames(d.bak), colnames(m.bak))

tb = tb[match(colnames(d.bak), tb[,1]), ]
scalematrix <- function(data) {
  cm <- rowMeans(data)
  csd <- sqrt((rowMeans(data*data) - cm^2) / (ncol(data) - 1) * ncol(data))
  (data - cm) / csd
}

corfunc <- function(m1,m2,type='concordant') {
  rowSums(scalematrix(m1) * scalematrix(m2))/(ncol(m1)-1)
}

sfunc <- function(tmpm) {
  cm <- rowMeans(tmpm) ## CpG site mean
  csd <- sqrt((rowMeans(tmpm*tmpm) - cm^2) / (ncol(tmpm) - 1) * ncol(tmpm)) ## CpG site sd
  names(head(sort(csd,decreasing = T),100000))  
}


for (disease in rev(unique(tb[,3]))){
  id = tb[tb[,3] == disease, 1]
  str(id)
  d = d.bak[, id]
  m = m.bak[, id]
  
  ## ======================================= ##
  ## use matched RNA-DNAm for model training ## 
  ## ======================================= ##
  set.seed(12345)
  trainid <- sample(colnames(m),ncol(m)*0.8)
  testid <- setdiff(colnames(m),trainid)
  saveRDS(trainid, paste0(rdir, 'trainid_', ifelse(disease, 'cancer', 'normal'), '.rds')) 
  saveRDS(testid, paste0(rdir, 'testid_', ifelse(disease, 'cancer', 'normal'), '.rds')) 
  
  sitevar = sapply(1:nrow(m), function(i){
    var(m[i, testid])
  })
  saveRDS(sitevar, paste0(rdir, 'sitevar_in_testid_', ifelse(disease, 'cancer', 'normal'), '.rds')) 
  
  sitevar = sapply(1:nrow(m), function(i){
    var(m[i, trainid])
  })
  saveRDS(sitevar, paste0(rdir, 'sitevar_in_trainid_', ifelse(disease, 'cancer', 'normal'), '.rds')) 

  hsdid <- sfunc(m[,testid]) ## select high sd CpG sites in testing side
  mod <- trainmodel(d[,trainid],m[,trainid])
  pred <- predict(d[,testid],mod)
  saveRDS(pred, paste0(rdir, 'predicted_DNAm_on_testset_', ifelse(disease,'cancer', 'normal'),'.rds'))
  
  true <- m[,testid]
  sampcv.all <- corfunc(pred,true) ## cross-sample PCC, only for highsd CpG
  saveRDS(sampcv.all, paste0(rdir, 'sampcvall_', ifelse(disease, 'cancer', 'normal'), '.rds'))
  
  sitecv <- corfunc(t(pred),t(true)) ## cross-CpG PCC
  sampcv <- sampcv.all[hsdid]## cross-sample PCC, only for highsd CpG
  saveRDS(sampcv, paste0(rdir, 'sampcv_', ifelse(disease, 'cancer', 'normal'), '.rds'))
  saveRDS(sitecv, paste0(rdir, 'sitecv_', ifelse(disease, 'cancer', 'normal'), '.rds'))
  
  summary(sampcv)
  summary(sitecv)
  
  mse <- sapply(1:nrow(pred), function(i){
    tmpv <- pred[i, ] - true[i, ]
    mean(tmpv * tmpv)
  })
  
  names(mse) <- rownames(pred)
  saveRDS(mse, paste0(rdir, 'mse_', ifelse(disease,'cancer', 'normal'),'.rds'))

  ## ======================================================== ##
  ## permute the RNA-DNAm relationship, and retrain the model ## 
  ## ======================================================== ##
  set.seed(12345)
  mod <- trainmodel(d[,trainid],m[,sample(trainid)]) 
  pred <- predict(d[,testid],mod)
  saveRDS(pred, paste0(rdir, 'predicted_DNAm_on_testset_in_permuted_model_', ifelse(disease,'cancer', 'normal'),'.rds'))
  
  true <- m[,testid]
  persampcvall <- corfunc(pred,true)
  saveRDS(persampcvall, paste0(rdir, 'persampcvall_', ifelse(disease, 'cancer', 'normal'), '.rds'))
  
  persampcv <- persampcvall[hsdid]
  persitecv <- corfunc(t(pred),t(true))
  saveRDS(persampcv, paste0(rdir, 'persampcv_', ifelse(disease, 'cancer', 'normal'), '.rds'))
  saveRDS(persitecv, paste0(rdir, 'persitecv_', ifelse(disease, 'cancer', 'normal'), '.rds'))
  summary(persampcv)
  summary(persitecv)

  permse <- sapply(1:nrow(pred), function(i){
    tmpv <- pred[i, ] - true[i, ]
    mean(tmpv * tmpv)
  })
  names(permse) <- rownames(pred)
  saveRDS(permse, paste0(rdir, 'permse_', ifelse(disease,'cancer', 'normal'),'.rds'))

    
  
  pd <- data.frame(cv=c(sampcv,sitecv,persampcv,persitecv),type=rep(c('across-sample','across-loci','across-sample\n(unmatched)','across-loci\n(unmatched)'),c(length(sampcv),length(sitecv),length(persampcv),length(persitecv))),stringsAsFactors = F)
  
  library(ggplot2)
  pdf(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/metpred/predict/plot/bulkcv/acrosssamplecor_', ifelse(disease, 'cancer', 'normal'), '.pdf'),width=2.7,height=1.8)
  print(ggplot(pd[pd$type%in%c('across-sample','across-sample\n(unmatched)'),],aes(cv,type,col=type, fill = type)) + geom_violin(alpha = 0.5, scale = 'width') + scale_color_manual(values=c('across-sample'='royalblue','across-sample\n(unmatched)'='grey'))  + scale_fill_manual(values=c('across-sample'='royalblue','across-sample\n(unmatched)'='grey')) + geom_boxplot(width=0.2,  outlier.size = 0.3) + coord_flip() + theme_classic() + theme(legend.position = 'none', axis.text.x = element_text(color = 'black')) + ylab('') + xlab('Correlation'))
  dev.off()
  
  pdf(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/metpred/predict/plot/bulkcv/acrosscitecor_', ifelse(disease, 'cancer', 'normal'), '.pdf'),width=2.7,height=1.8)
  print(ggplot(pd[pd$type%in%c('across-loci','across-loci\n(unmatched)'),],aes(cv,type,col=type, fill = type)) + geom_violin(alpha = 0.5, scale = 'width') + scale_color_manual(values=c('across CpG sites'='orange','across-loci\n(unmatched)'='grey')) + geom_boxplot(width=0.2, outlier.size = 0.3) + coord_flip() + theme_classic() + theme(legend.position = 'none', axis.text.x = element_text(color = 'black')) + ylab('') + xlab('Correlation'))
  dev.off()
  
  id <- which.max(sampcv)
  pdf(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/metpred/predict/plot/bulkcv/acrosssampleexp_', ifelse(disease, 'cancer', 'normal'), '.pdf'), width=2.4,height=2.1)
  print(ggplot(data=data.frame(x=pred[id,],y=true[id,]),aes(x=x,y=y)) + geom_point(col='royalblue') + theme_classic() + xlab('Predicted methylation') + ylab('True methylation') + theme(legend.position = 'none') + geom_abline(slope=1,intercept = 0,col='red',linetype=2))
  dev.off()
  
  id <- which.max(sitecv)
  samppid <- sample(1:nrow(pred),200000)
  pdf(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/metpred/predict/plot/bulkcv/acrossciteexp_', ifelse(disease, 'cancer', 'normal'), '.pdf'),width=2.4,height=2.1)
  print(ggplot(data=data.frame(x=pred[samppid,id],y=true[samppid,id]),aes(x=x,y=y)) + geom_point(alpha=0.1,size=0.1,col='orange') + theme_classic() + xlab('Predicted methylation') + ylab('True methylation') + theme(legend.position = 'none') + geom_abline(slope=1,intercept = 0,col='red',linetype=2))
  dev.off()
  
}



