source('/home/whou10/data/whou/metpred/software/trainpredict.R')
source('/home/whou10/scratch16/whou10/resource/startup.R')
ddir <- '/home/whou10/data/whou/metpred/data/pbmc/imp/res/'
jid <- as.numeric(commandArgs(trailingOnly = T))
print(jid)
me = readRDS(paste0(ddir, 'wgbs_train.rds'))
rgroup <- as.numeric(cut(1:nrow(me),20))
me = me[rgroup==jid,]
me = me[!sub(':.*','',rownames(me)) %in% c('chrM','chrY'),]
r.tr = readRDS(paste0(ddir, 'rna_train.rds'))
r.te = readRDS(paste0(ddir, 'rna_test_sc.rds'))
ct = sub(':.*', '', colnames(r.te))

r.te = r.te[, ct %in% c('b_cells', 'cd14_monocytes', 'cd56_nk', 
                        'cytotoxic_t', 'memory_t', 'naive_t',
                        'regulatory_t') ]


ct = sub(':.*', '', colnames(r.te))
id <- unlist(sapply(unique(ct), function(uct){
  which(ct == uct)[1:300]
}, simplify = F))
r.te = r.te[,id]

# new program
source('/home/whou10/data/whou/metpred/software/nearestgenepredict.R')




nearestgenepredict <- function(trainexpr,testexpr,trainmeth,species='human') {
  library(data.table)
  library(fastcluster)
  library(reshape2)
  suppressMessages(library(GenomicRanges))
  
  logistic <- function(data) {
    1/(1+exp(-data))
  }
  
  rowsds <- function(data,cm) {
    sqrt((rowMeans(data*data) - cm^2) / (ncol(data) - 1) * ncol(data))  
  }
  
  trainmeth <- trainmeth[,colnames(trainexpr)]
  trainmeth[trainmeth==0] <- min(trainmeth[trainmeth>0])
  trainmeth[trainmeth==1] <- max(trainmeth[trainmeth<1])
  
  trainmeth <- log(trainmeth/(1-trainmeth))
  
  int <- intersect(rownames(trainexpr),rownames(testexpr))
  trainexpr <- trainexpr[int,]
  testexpr <- testexpr[int,]
  colnames(trainmeth) <- colnames(trainexpr) <- paste0('train_',1:ncol(trainexpr))
  oritestname <- colnames(testexpr)
  names(oritestname) <- colnames(testexpr) <- paste0('test_',1:ncol(testexpr))
  
  expr <- cbind(trainexpr,testexpr)
  dn <- dimnames(expr)
  qnem <- rowMeans(apply(expr,2,function(i) sort(i)))
  expr <- apply(expr,2,function(i) qnem[frank(i)])
  dimnames(expr) <- dn
  
  trainexpr <- expr[,paste0('train_',1:ncol(trainexpr))]
  testexpr <- expr[,paste0('test_',1:ncol(testexpr))]
  
  trainexpr <- trainexpr[rowMeans(trainexpr > median(trainexpr[trainexpr > 0])) >= 0.1,]
  m <- rowMeans(trainexpr)
  s <- rowsds(trainexpr,m)
  trainexpr <- trainexpr[s/m > 0.1,]
  
  if (!grepl('-',rownames(me)[1])) {
    megrchr <- sub(':.*','',rownames(me))
    megrpos <- as.numeric(sub('.*:','',rownames(me)))
    megr <- GRanges(seqnames=megrchr,IRanges(start=megrpos,end=megrpos))
  } else {
    megrchr <- sub(':.*','',rownames(me))
    megrstart <- as.numeric(sub('-.*','',sub('.*:','',rownames(me))))
    megrend <- as.numeric(sub('.*-','',sub('.*:','',rownames(me))))
    megr <- GRanges(seqnames=megrchr,IRanges(start=megrstart,end=megrend))
  }
  print(megr)
  
  if (species=='human') {
    gtf <- data.table::fread('/home/whou10/data/whou/resource/gtf/grch38.gtf',data.table=F)  
  } else if (species=='mouse') {
    gtf <- data.table::fread('/home/whou10/data/whou/resource/gtf/grcm38.gtf',data.table=F)  
  }
  gtf <- gtf[gtf[,3]=='gene',]
  gn <- sub('\".*','',sub('.*gene_name "','',gtf[,9]))
  trainexpr <- trainexpr[rownames(trainexpr) %in% gn,]
  
  tss <- GRanges(seqnames=gtf[,1],IRanges(start=gtf[,4],end=gtf[,5]),strand=gtf[,7])
  tss <- promoters(tss,1,0)
  names(tss) <- gn
  tss <- tss[rownames(trainexpr)]
  
  ne <- names(tss)[nearest(megr,tss)]
  ne <- split(1:length(ne),ne)
  pred <- sapply(names(ne),function(g) {
    x <- cbind(1,trainexpr[g,])
    tx <- cbind(1,testexpr[g,])
    trainmeth[ne[[g]],,drop=F]  %*% (x %*% chol2inv(chol(crossprod(x))) %*% t(tx))
  },simplify = F)
  
  pred <- do.call(rbind,pred)
  pred <- logistic(pred)
  colnames(pred) <- oritestname[colnames(pred)]
  pred[rownames(trainmeth),]
}




rownames(me) <- sub('_',':',rownames(me))
pred <- nearestgenepredict(trainexpr=r.tr,testexpr=r.te,trainmeth=me   )
print('done')
saveRDS(pred,file=paste0('/home/whou10/data/whou/metpred/evaluate/pbmc/predsc/imp/competing/res/neargene_',jid,'.rds'))


# samplecv <- sapply(1:ncol(pred),function(i) cor(pred[,i],testme[,i]))
# sitecv <- sapply(1:nrow(pred),function(i) cor(pred[i,],testme[i,]))
# testsd <- apply(testme,1,sd)






