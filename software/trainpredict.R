trainpredict <- function(trainexpr,testexpr,trainmeth) {
  library(data.table)
  library(fastcluster)
  library(reshape2)
  
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
  stdtrainexpr <- (trainexpr-m)/s
  stdtrainexpr <- stdtrainexpr[s/m > 0.1,]
  
  clunumlist <- c(1000,2500,5000)
  clunumlist <- clunumlist[clunumlist < nrow(stdtrainexpr)]
  cluexpr <- sapply(clunumlist,function(cn) {
    set.seed(1)
    gclu <- kmeans(stdtrainexpr,cn,iter.max=10000)$cluster
    cluexpr <- rowsum(expr[names(gclu),],gclu)
    cluexpr/as.vector(table(gclu)[rownames(cluexpr)])
  },simplify = F)
  names(cluexpr) <- clunumlist
  
  #intrainexpr=cvtrainexpr;intestexpr=cvtestexpr;intrainmeth=cvtrainmeth;clunum=5;lambdalist=0
  
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
  
  ### cross validation
  set.seed(12345)
  spid <- split(sample(colnames(trainexpr)),ceiling(seq_along(colnames(trainexpr))/(ncol(trainexpr)/5)))
  perf <- NULL
  for (i in 1:length(spid)) {
    print(paste0('CV fold ',i))
    cvtestid <- spid[[i]]
    cvtrainid <- setdiff(colnames(trainexpr),cvtestid)
    cvtrainmeth <- trainmeth[,cvtrainid]
    cvtestmeth <- logistic(trainmeth[,cvtestid])
    
    lambdalist <- c(10^c(-1,0,1,2))
    for (clunum in clunumlist) {
      pred <- inpredfunc(cluexpr[[as.character(clunum)]][,cvtrainid],cluexpr[[as.character(clunum)]][,cvtestid],cvtrainmeth,lambdalist)
      mse <- sapply(pred,function(i) rowMeans((i-cvtestmeth)^2))
      colnames(mse) <- lambdalist
      mse <- melt(mse)
      colnames(mse) <- c('cpg','lambda','mse')
      perf <- rbind(perf,data.frame(clunum=clunum,mse))
    }
  }
  perf$combpara <- paste0(perf$clunum,'_',perf$lambda)
  perf <- aggregate(perf$mse,list(perf$combpara,perf$cpg),mean)
  perf <- perf[order(perf$x),]
  perf <- perf[!duplicated(perf[,2]),]
  
  print(table(perf[,1]))
  pred <- sapply(unique(perf[,1]),function(uset) {
    id <- which(perf[,1]==uset)
    clunum=sub('_.*','',uset)
    inpredfunc(cluexpr[[clunum]][,colnames(trainexpr)],cluexpr[[clunum]][,colnames(testexpr)],trainmeth[id,,drop=F],lambda=as.numeric(sub('.*_','',uset)))
  },simplify = F)
  pred <- do.call(rbind,pred)
  
  colnames(pred) <- oritestname[colnames(pred)]
  pred[rownames(trainmeth),]
}


