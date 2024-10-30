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
  trainmeth[trainmeth==0] <- 1e-3
  trainmeth[trainmeth==1] <- 1-1e-3
  
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
  
  if (!grepl('-',rownames(trainmeth)[1])) {
    megrchr <- sub(':.*','',rownames(trainmeth))
    megrpos <- as.numeric(sub('.*:','',rownames(trainmeth)))
    megr <- GRanges(seqnames=megrchr,IRanges(start=megrpos,end=megrpos))
  } else {
    megrchr <- sub(':.*','',rownames(trainmeth))
    megrstart <- as.numeric(sub('-.*','',sub('.*:','',rownames(trainmeth))))
    megrend <- as.numeric(sub('.*-','',sub('.*:','',rownames(trainmeth))))
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
  pred[intersect(rownames(trainmeth), rownames(pred)),]
}

