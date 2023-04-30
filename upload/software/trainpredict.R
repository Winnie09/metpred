trainpredict <- function(trainexpr,testexpr,trainmeth) {
  library(data.table)
  library(fastcluster)
  library(reshape2)
  
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
  
  #intrainexpr=cvtrainexpr;intestexpr=cvtestexpr;intrainmeth=cvtrainmeth;clunum=5;lambdalist=0
  
  inpredfunc <- function(intrainexpr,intestexpr,intrainmeth,clunum,lambdalist) {
    intrainexpr <- intrainexpr[rowMeans(intrainexpr > median(intrainexpr[intrainexpr > 0])) >= 0.1,]
    m <- rowMeans(intrainexpr)
    s <- sqrt((rowMeans(intrainexpr*intrainexpr) - m^2) * ncol(intrainexpr)/(ncol(intrainexpr)-1))
    stdintrainexpr <- (intrainexpr-m)/s
    stdintrainexpr <- stdintrainexpr[s/m > 1,]
    
    set.seed(1)
    gclu <- kmeans(stdintrainexpr,clunum,iter.max=10000)$cluster
    
    cluintrainexpr <- rowsum(intrainexpr[names(gclu),],gclu)
    cluintrainexpr <- cluintrainexpr/as.vector(table(gclu)[rownames(cluintrainexpr)])
    cluintestexpr <- rowsum(intestexpr[names(gclu),],gclu)
    cluintestexpr <- cluintestexpr/as.vector(table(gclu)[rownames(cluintestexpr)])
    
    cluintrainexprmean <- rowMeans(cluintrainexpr)
    cluintrainexprsd <- apply(cluintrainexpr,1,sd)
    
    cluintrainexprscale <- t((cluintrainexpr-cluintrainexprmean)/cluintrainexprsd)
    cluintestexprscale <- t((cluintestexpr-cluintrainexprmean)/cluintrainexprsd)
    crossx <- crossprod(cluintrainexprscale)
    
    intrainmethmean <- rowMeans(intrainmeth)
    intrainmethsd <- apply(intrainmeth,1,sd)
    intrainmethscale <- (intrainmeth-intrainmethmean)/intrainmethsd
    
    out <- sapply(lambdalist,function(lambda) {
      diagm <- diag(rep(lambda,ncol(cluintrainexprscale)))*nrow(cluintrainexprscale)
      res <- intrainmethscale %*% (cluintrainexprscale %*% chol2inv(chol(crossx+diagm)) %*% t(cluintestexprscale)) * intrainmethsd + intrainmethmean
      if (sum(intrainmethsd==0) > 0)
        res[intrainmethsd==0,] <- intrainmethmean[intrainmethsd==0]
      res
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
    cvtrainexpr <- trainexpr[,cvtrainid]
    cvtestexpr <- trainexpr[,cvtestid]
    cvtrainmeth <- trainmeth[,cvtrainid]
    cvtestmeth <- trainmeth[,cvtestid]
    
    lambdalist <- c(10^c(-1,0,1,2,3))
    for (clunum in c(1000,2500,5000)) {
      pred <- inpredfunc(cvtrainexpr,cvtestexpr,cvtrainmeth,clunum,lambdalist)
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
    inpredfunc(trainexpr,testexpr,trainmeth[id,,drop=F],clunum=as.numeric(sub('_.*','',uset)),lambda=as.numeric(sub('.*_','',uset)))  
  },simplify = F)
  pred <- do.call(rbind,pred)
  
  pred <- exp(pred)
  pred <- pred/(1+pred)
  colnames(pred) <- oritestname[colnames(pred)]
  pred[rownames(trainmeth),]
}


