## these codes are used to debug
## why some cpg was predicted as very low, but the gs is very high

rm(list=ls())
## prepare evaluation data: goldstandard, and predicted sc pooled pb
## for all cpg prediction
source('/home/whou10/scratch16/whou10/resource/myfunc/01_function.R')
ddir <- '/home/whou10/data/whou/metpred/data/pbmc/imp/res/'
rdir = '/home/whou10/data/whou/metpred/evaluate/pbmc/predsc/imp/res/eval/'
w.te = readRDS(paste0(ddir, 'wgbs_test.rds'))
pred.agg = readRDS('/home/whou10/data/whou/metpred/evaluate/pbmc/predsc/imp/res/pb/pred_allcpg.rds')
w.te2 = w.te[rownames(pred.agg),]
w.te.agg = aggregatefunc2(d = w.te2, 
                          by = c('cd14_monocytes', 'b_cells', rep('t_cells', 3)), 
                          fun = 'mean')
pred.agg = pred.agg[, colnames(w.te.agg)]
# saveRDS(pred.agg, paste0(rdir, 'eval_pred_ctmean.rds'))
# saveRDS(w.te.agg, paste0(rdir, 'eval_goldstandard_ctmean.rds'))

ng = readRDS('/home/whou10/data/whou/metpred/evaluate/pbmc/predsc/imp/competing/combine/neargene.rds')
str(ng)

pdir = '/home/whou10/data/whou/metpred/evaluate/pbmc/predsc/imp/plot/'
# png(paste0(pdir, 'compare_pred_and_gs_smoothscatter.png'), width = 2e3, height = 2e3, res = 200)
# par(mfrow=c(3,2))
# mdf = data.frame(v1 = c(1,2,3,4,5), v2=c(2,1,4,4,4))
# for (i in 1:nrow(mdf)){
#   print(i)
#   smoothScatter(w.te2[,mdf[i,1]], 
#                 pred.agg[,mdf[i,2]],
#                 xlab=colnames(w.te2)[mdf[i,1]], 
#                 ylab=colnames(pred.agg)[mdf[i,2]], 
#                 main = round(cor(w.te2[,mdf[i,1]], pred.agg[,mdf[i,2]]), 3))
#   
# }
# dev.off()


int = intersect(rownames(ng), rownames(w.te.agg))
png(paste0(pdir, 'compare_neargene_and_gs_smoothscatter.png'), width = 2e3, height = 2e3, res = 200)
par(mfrow=c(2,2))
mdf = data.frame(v1 = c(1,2,3), v2=c(2,1,4))
for (i in 1:nrow(mdf)){
  print(i)
  smoothScatter(w.te.agg[int, mdf[i,1]],
                ng[int, mdf[i,2]],
                xlab=colnames(w.te.agg)[mdf[i,1]],
                ylab=colnames(ng)[mdf[i,2]],
                main = round(cor(w.te.agg[int,mdf[i,1]], ng[int,mdf[i,2]]), 3))

}
dev.off()

## check why: very high in gs, very low in pred
id <- which(pred.agg[,1] < 0.05 & w.te.agg[,1] > 0.95)
id2 <- which(pred.agg[,2] < 0.05 & w.te.agg[,2] > 0.95)
id3 <- which(pred.agg[,3] < 0.05 & w.te.agg[,3] > 0.95)
idint = intersect(intersect(id, id2), id3)
str(idint)
names(idint) = rownames(pred.agg)[idint]

str(pred.agg)
str(w.te.agg)
png(paste0(pdir, 'compare_pred_and_gs_smoothscatter_cpgidint.png'), width = 2e3, height = 2e3, res = 200)
par(mfrow=c(2,2))
mdf = data.frame(v1 = c(1,2,3), v2=c(1,2,3))
for (i in 1:nrow(mdf)){
  print(i)
  smoothScatter(w.te.agg[idint, mdf[i,1]],
                pred.agg[idint, mdf[i,2]],
                xlab=colnames(w.te.agg)[mdf[i,1]],
                ylab=colnames(pred.agg)[mdf[i,2]],
                main = round(cor(w.te.agg[idint,mdf[i,1]], pred.agg[idint,mdf[i,2]]), 3))
  
}
dev.off()


# > str(id)
# Named int [1:306256] 539697 539698 539699 539700 539701 539702 539703 539704 539705 539706 ...
# - attr(*, "names")= chr [1:306256] "chr1_85707895" "chr1_85707922" "chr1_85707945" "chr1_85707951" ...
pheatmap(apply(w.te.agg[id,],2,summary), scale = 'none')
# cd14_monocytes       b_cells      t_cells
# Min.         0.9500002 5.026057e-124 8.714054e-06
# 1st Qu.      0.9546530  8.621937e-01 8.809705e-01
# Median       0.9593355  9.099717e-01 9.241050e-01
# Mean         0.9605258  8.888304e-01 9.049166e-01
# 3rd Qu.      0.9649428  9.449102e-01 9.532470e-01
# Max.         1.0000000  1.000000e+00 1.000000e+00
pheatmap(apply(pred.agg[id,],2,summary), scale = 'none')
# cd14_monocytes       b_cells       t_cells
# Min.     4.832109e-156 5.093174e-146 1.449595e-154
# 1st Qu.   1.320688e-02  1.003149e-02  1.690018e-03
# Median    2.138835e-02  1.664814e-02  2.828271e-03
# Mean      2.285463e-02  1.862363e-02  3.323421e-03
# 3rd Qu.   3.178620e-02  2.554848e-02  4.321837e-03
# Max.      4.999968e-02  9.979958e-01  6.281201e-01
w.tr = readRDS(paste0(ddir, 'wgbs_train.rds'))
w.tr = w.tr[rownames(w.te.agg),]
str(w.tr)
# does these cpg has special patterns in training DNAm?
library(pheatmap)
m = apply(w.tr[id,],2,summary)
pheatmap(m, scale='none')
rm_te = rowMeans(w.te.agg[id,])
rm_tr = rowMeans(w.tr[id,])
smoothScatter(rm_te, rm_tr) ## quite consistent 
## cpg: high in test, low in training 
id2 = which(rm_te > 0.95 & rm_tr < 0.05)
# > str(id2)
# Named int(0) 
# - attr(*, "names")= chr(0) 
## they are well predicted

## check in pred and gs
png(paste0(pdir, 'compare_pred_and_gs_smoothscatter_cpgid.png'), width = 2e3, height = 2e3, res = 200)
par(mfrow=c(3,2))
mdf = data.frame(v1 = c(1,2,3,4,5), v2=c(1,2,3,3,3))
for (i in 1:nrow(mdf)){
  print(i)
  smoothScatter(w.te2[names(id),mdf[i,1]],
                pred.agg[names(id),mdf[i,2]],
                xlab=colnames(w.te2)[mdf[i,1]],
                ylab=colnames(pred.agg)[mdf[i,2]],
                main = round(cor(w.te2[names(id),mdf[i,1]], pred.agg[names(id),mdf[i,2]]), 3))

}
dev.off()

## compare training and testing. on terribly-predicted cpg sites 
png(paste0(pdir, 'compare_tr_and_te_smoothscatter_cpgid.png'), width = 2e3, height = 2e3, res = 200)
par(mfrow=c(3,3))
mdf = rbind(data.frame(v1 = grep('monocyte', colnames(w.tr)),
           v2 = 1, stringsAsFactors = F),
data.frame(v1 = grep('B_cell', colnames(w.tr)),
           v2 = 2, stringsAsFactors = F),
data.frame(v1 = grep('T_cell', colnames(w.tr)),
           v2 = 3, stringsAsFactors = F))
for (i in 1:nrow(mdf)){
  print(i)
  smoothScatter(w.tr[names(id),mdf[i,1]],
                w.te2[names(id),mdf[i,2]],
                xlab=colnames(w.tr)[mdf[i,1]],
                ylab=colnames(w.te2)[mdf[i,2]],
                main = round(cor(w.tr[names(id),mdf[i,1]], w.te2[names(id),mdf[i,2]]), 3),
                ylim=c(0,1),
                xlim=c(0,1))
  
}
dev.off()


## check if neargene also has terrible prediction on these cpg id
png(paste0(pdir, 'compare_neargene_and_gs_smoothscatter_cpgid.png'), width = 2e3, height = 2e3, res = 200)
par(mfrow=c(2,2))
mdf = data.frame(v1 = c(1,2,3), v2=c(2,1,4))
v = names(id)[names(id) %in% int]
str(v)
for (i in 1:nrow(mdf)){
  print(i)
  smoothScatter(w.te.agg[v, mdf[i,1]],
                ng[v, mdf[i,2]],
                xlab=colnames(w.te.agg)[mdf[i,1]],
                ylab=colnames(ng)[mdf[i,2]],
                main = round(cor(w.te.agg[v,mdf[i,1]], ng[v,mdf[i,2]]), 3))
  
}
dev.off()



## ===============================================
## check through the prediction codes
## ===============================================
set.seed(1)
cpg = names(idint)[sample(1:length(idint), 100)]



rs2 = readRDS(paste0(ddir, 'wgbs_train_nonNA_CpG_rowsds_cut.rds'))
str(rs2)
table(rs2)
# rs2
# (0,0.05] (0.05,0.1] (0.1,0.15] (0.15,0.2] (0.2,0.25] (0.25,0.3] 
# 12290958   12632061    4832776     592779     203745      71496 
# (0.3,0.35] (0.35,0.4] (0.4,0.45] (0.45,0.5] (0.5,0.55] (0.55,0.6] 
# 29446      14198       4619        904          0          0 
# (0.6,0.65] (0.65,0.7] (0.7,0.75] (0.75,0.8] (0.8,0.85] (0.85,0.9] 
# 0          0          0          0          0          0 
# (0.9,0.95]   (0.95,1] 
# 0          0 
table(rs2[names(idint)])
# (0,0.05] (0.05,0.1] (0.1,0.15] (0.15,0.2] (0.2,0.25] (0.25,0.3] 
# 189377      90318      20621       2868       1047        433 
# (0.3,0.35] (0.35,0.4] (0.4,0.45] (0.45,0.5] (0.5,0.55] (0.55,0.6] 
# 590        552        384         66          0          0 
# (0.6,0.65] (0.65,0.7] (0.7,0.75] (0.75,0.8] (0.8,0.85] (0.85,0.9] 
# 0          0          0          0          0          0 
# (0.9,0.95]   (0.95,1] 
# 0          0 
# w.tr = w.tr[names(rs2)[rs2 == levels(rs2)[cpggroup]], , drop = F]
table(rs2[cpg])
trainmeth = w.tr[names(idint), , drop = F]
# > str(trainmeth)
# num [1:100, 1:66] 0.948 0.894 0.967 0.973 0.972 ...
# - attr(*, "dimnames")=List of 2
# ..$ : chr [1:100] "chr1_202769009" "chr1_151315880" "chr1_150958056" "chr1_161913292" ...
# ..$ : chr [1:66] "blueprint_venous_blood-C000S5-CD14-positive_CD16-negative_classical_monocyte" "blueprint_venous_blood-C000S5-mature_neutrophil" "blueprint_venous_blood-C0010K-CD14-positive_CD16-negative_classical_monocyte" "blueprint_venous_blood-C0010K-mature_neutrophil" ...
if (nrow(w.tr) < 1){
  stop('No CpG in this rowsds group!')
}

r.tr = readRDS(paste0(ddir, 'rna_train.rds'))
str(r.tr)
# num [1:58434, 1:66] 0 0 2.087 0.566 0.506 ...
# - attr(*, "dimnames")=List of 2
# ..$ : chr [1:58434] "TSPAN6" "TNMD" "DPM1" "SCYL3" ...
# ..$ : chr [1:66] "blueprint_venous_blood-C000S5-CD14-positive_CD16-negative_classical_monocyte" "blueprint_venous_blood-C000S5-mature_neutrophil" "blueprint_venous_blood-C0010K-CD14-positive_CD16-negative_classical_monocyte" "blueprint_venous_blood-C0010K-mature_neutrophil" ...
r.te.sub = readRDS(paste0(ddir, 'rna_test_sc_sub300.rds'))
# > str(r.te.sub)
# num [1:8209, 1:2100] 0.0963 0.3414 0.0868 0.0426 0.1903 ...
# - attr(*, "dimnames")=List of 2
# ..$ : chr [1:8209] "NOC2L" "ISG15" "TNFRSF18" "TNFRSF4" ...
# ..$ : chr [1:2100] "b_cells:AAACATACAATGCC-1" "b_cells:AAACATACACGCAT-1" "b_cells:AAACATACGAATAG-1" "b_cells:AAACATACGTGTCA-1" ...
# ct = sub(':.*', '', colnames(r.te.sub))
# table(ct)
# b_cells cd14_monocytes        cd56_nk    cytotoxic_t       memory_t 
# 300            300            300            300            300 
# naive_t   regulatory_t 
# 300            300 
trainexpr=r.tr
testexpr=r.te.sub
trainmeth=trainmeth
clunumlist = c(1e3)
lambdalist = c(10e-1)


trainpredict <- function(trainexpr,testexpr,trainmeth,clunumlist = c(1000,2500,5000), lambdalist = c(10^c(-1,0,1))) {
  ## if the lengths of clunumlist and lambdalist are both 1, no cross validation
  ## if lengths of clunumlist and lambdalist are not both 1 (at least one length > 1), then conduct cross validation to optimize the values. 
  library(data.table)
  library(fastcluster)
  library(reshape2)
  
  logistic <- function(data) {
    1/(1+exp(-data))
  }
  
  rowsds <- function(data,cm) {
    sqrt((rowMeans(data*data) - cm^2) / (ncol(data) - 1) * ncol(data))  
  }
  
  
  
  complete_cases_rowwise <- function(mat) {
    apply(mat, 1, function(row) all(!is.na(row)))
  }
  
  # Apply the function to your data
  complete_cases <- complete_cases_rowwise(trainmeth)
  
  # Subset the data using the complete cases index
  trainmeth_complete <- trainmeth[complete_cases, , drop = FALSE]
  
  # Check the dimensions to ensure it worked
  trainmeth <- trainmeth_complete
  print(paste0('Reserved ', nrow(trainmeth), ' CpGs.'))
  
  # ### note: long vectors not supported yet: complete_cases.c:192
  # if (sum(complete.cases(trainmeth)) != nrow(trainmeth)) {
  #   print('Training methylation contains NA. Removing rows with NAs.')
  #   trainmeth = trainmeth[complete.cases(trainmeth), , drop = F]
  #   print(paste0('Reserved ', nrow(trainmeth), ' CpGs.'))
  # }
  
  trainmeth <- trainmeth[,colnames(trainexpr)]
  # trainmeth[trainmeth==0] <- min(trainmeth[trainmeth>0]) ## 0.4098026
  # trainmeth[trainmeth==1] <- max(trainmeth[trainmeth<1]) ## 0.9939805
  trainmeth[trainmeth==0] <- 1e-3
  trainmeth[trainmeth==1] <- 1-1e-3
  
  ## training mean, testing mean
  
  ## logit transformation: from (0,1) to (-inf, inf). Need this for Ridge regression.
  trainmeth <- log(trainmeth/(1-trainmeth))
  
  ## use intereception of genes in train and test expr, mark train and test in methylation column names
  int <- intersect(rownames(trainexpr),rownames(testexpr))
  print(paste0('Number of overlapping genes in train and test expr: ', length(int)))
  trainexpr <- trainexpr[int,,drop=FALSE]
  testexpr <- testexpr[int,,drop=FALSE]
  colnames(trainmeth) <- colnames(trainexpr) <- paste0('train_',1:ncol(trainexpr))
  oritestname <- colnames(testexpr)
  names(oritestname) <- colnames(testexpr) <- paste0('test_',1:ncol(testexpr))
  
  ## quantile normalize train and test expr, mark down train and test in expr column names
  expr <- cbind(trainexpr,testexpr)
  dn <- dimnames(expr)
  qnem <- rowMeans(apply(expr,2,function(i) sort(i)))
  expr <- apply(expr,2,function(i) qnem[frank(i)])
  dimnames(expr) <- dn
  
  trainexpr <- expr[,paste0('train_',1:ncol(trainexpr)),drop=FALSE]
  testexpr <- expr[,paste0('test_',1:ncol(testexpr)),drop=FALSE]
  
  trainexpr <- trainexpr[rowMeans(trainexpr > median(trainexpr[trainexpr > 0])) >= 0.1,,drop=FALSE]
  m <- rowMeans(trainexpr)
  s <- rowsds(trainexpr,m)
  stdtrainexpr <- (trainexpr-m)/s
  stdtrainexpr <- stdtrainexpr[s/m > 0.1,]
  
  ## gene cluster number
  clunumlist <- clunumlist[clunumlist < nrow(stdtrainexpr)]
  cluexpr <- sapply(clunumlist,function(cn) {
    set.seed(1)
    gclu <- kmeans(stdtrainexpr,cn,iter.max=10000)$cluster
    cluexpr <- rowsum(expr[names(gclu),],gclu)
    cluexpr/as.vector(table(gclu)[rownames(cluexpr)])
  },simplify = F)
  names(cluexpr) <- clunumlist
  
  #cluintrainexpr=cluexpr[[clunum]][,colnames(trainexpr)];cluintestexpr=cluexpr[[clunum]][,colnames(testexpr)];intrainmeth=trainmeth;lambdalist=lambdalist
  
  inpredfunc <- function(cluintrainexpr,cluintestexpr,intrainmeth,lambdalist) {
    print(1)
    cluintrainexprmean <- rowMeans(cluintrainexpr)
    cluintrainexprsd <- rowsds(cluintrainexpr,rowMeans(cluintrainexpr))
    
    cluintrainexprscale <- t((cluintrainexpr-cluintrainexprmean)/cluintrainexprsd)
    cluintestexprscale <- t((cluintestexpr-cluintrainexprmean)/cluintrainexprsd)
    crossx <- crossprod(cluintrainexprscale)
    
    intrainmethmean <- rowMeans(intrainmeth)
    intrainmethsd <- rowsds(intrainmeth,intrainmethmean)
    intrainmethscale <- (intrainmeth-intrainmethmean)/intrainmethsd
    
    out <- sapply(lambdalist,function(lambda) {
      diagm <- diag(rep(lambda,ncol(cluintrainexprscale)))*nrow(cluintrainexprscale)
      res <- intrainmethscale %*% (cluintrainexprscale %*% chol2inv(chol(crossx+diagm)) %*% t(cluintestexprscale)) * intrainmethsd + intrainmethmean
      if (sum(intrainmethsd==0) > 0)
        res[intrainmethsd==0,] <- intrainmethmean[intrainmethsd==0]
      logistic(res)
    },simplify = F)  
    if (length(out)==1) {
      out[[1]]
    } else {
      out
    }
  }
  
  if (!length(clunumlist)==1 | !length(lambdalist)==1){
    print('Start cross validation...')
    set.seed(12345)
    spid <- split(sample(colnames(trainexpr)),ceiling(seq_along(colnames(trainexpr))/(ncol(trainexpr)/5)))
    perf <- NULL
    for (i in 1:length(spid)) {
      print(paste0('CV fold ',i))
      cvtestid <- spid[[i]]
      cvtrainid <- setdiff(colnames(trainexpr),cvtestid)
      cvtrainmeth <- trainmeth[,cvtrainid]
      cvtestmeth <- logistic(trainmeth[,cvtestid])
      
      for (clunum in clunumlist) { #####
        print(paste0('clunum ', clunum))	    
        pred <- inpredfunc(cluexpr[[as.character(clunum)]][,cvtrainid],cluexpr[[as.character(clunum)]][,cvtestid],cvtrainmeth,lambdalist)
        print('0...')
        print(str(pred))
        mse <- sapply(pred,function(i) rowMeans((i-cvtestmeth)^2))
        colnames(mse) <- lambdalist
        mse <- melt(mse)
        colnames(mse) <- c('cpg','lambda','mse')
        perf <- rbind(perf,data.frame(clunum=clunum,mse))
      }
    }
    perf$combpara <- paste0(perf$clunum,'_',perf$lambda)
    perf <- aggregate(perf$mse,list(perf$combpara,perf$cpg),mean)
    perf <- perf[order(perf$x),,drop=FALSE]
    print('1...')
    print(str(perf))
    
    perf <- perf[!duplicated(perf[,2]),,drop=FALSE]
    print('2...')
    print(str(perf)) 
    
    print('information...')
    print(table(perf[,1]))
    print(str(trainexpr))
    print(str(testexpr))
    pred <- sapply(unique(perf[,1]),function(uset) {
      id <- which(perf[,1]==uset)
      clunum=sub('_.*','',uset)
      inpredfunc(cluexpr[[clunum]][,colnames(trainexpr)],cluexpr[[clunum]][,colnames(testexpr)],trainmeth[id,,drop=F],lambda=as.numeric(sub('.*_','',uset)))
    },simplify = F)
    pred <- do.call(rbind,pred)
    print('End cross-validation.')
  } else { #######
    clunum = as.character(as.vector(clunumlist))
    pred = inpredfunc(cluexpr[[clunum]][,colnames(trainexpr)],cluexpr[[clunum]][,colnames(testexpr)],trainmeth,lambdalist)
  } 
  print(str(pred))
  colnames(pred) <- oritestname[colnames(pred)]
  #pred[rownames(trainmeth),]
  pred
}


## read in a prediction object that contains some of these terrible cpg
m = readRDS('/home/whou10/data/whou/metpred/evaluate/pbmc/predsc/imp/res/pred/allcpg/pred_cpggroup3.rds')
str(m)
int2 = intersect(rownames(m), rownames(pred))
v1 = m[int2,]
v2 = pred[int2,]
plot(v1, v2, xlab='predsave',ylab='preddebug')
st


str(r.tr)
str(trainmeth)
prednew = trainpredict(trainexpr=r.tr,
                       testexpr=r.te.sub,
                       trainmeth=w.tr[names(idint), , drop = F],
                       clunumlist = c(1e3),
                       lambdalist = c(10e-1))
saveRDS(list(trainexpr, testexpr, trainmeth, clunumlist, lambdalist), 
        '/home/whou10/data/whou/metpred/evaluate/pbmc/predsc/imp/res/pred/prednew_data.rds')
saveRDS(prednew, '/home/whou10/data/whou/metpred/evaluate/pbmc/predsc/imp/res/pred/prednew.rds')


summary(colMeans(prednew))
summary(colMeans(prednew[cpg, ]))
summary(colMeans(pred))
summary(colMeans(pred[cpg, ]))

summary(colMeans(prednew[rownames(pred), ]))



## ===============
## for J
## ===============
## new object 20241015: read in a prediction object that contains some of these terrible cpg
m = readRDS('/home/whou10/data/whou/metpred/evaluate/pbmc/predsc/imp/res/pred/allcpg/pred_cpggroup3.rds')

## old object 20240407: read in a prediction object that contains some of these terrible cpg
m2 = readRDS('/home/whou10/data/whou/metpred/evaluate/pbmc/predsc/imp/res/pred.bak/allcpg/pred_cpggroup3.rds')

## new object 20241014: predict on terrible cpg (low in prediction but high in goldstandard)
prednew = readRDS('/home/whou10/data/whou/metpred/evaluate/pbmc/predsc/imp/res/pred/prednew.rds')



int3 = intersect(rownames(m), rownames(prednew))
summary(colMeans(m[int3,]))
summary(colMeans(prednew[int3,]))
