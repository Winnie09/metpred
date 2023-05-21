expr = readRDS('/home/whou10/data/whou/metpred/software_test/rna_encode_sub_20.rds')
expr = expr[1:2e3,]
meth = readRDS('/home/whou10/data/whou/metpred/software_test/wgbs_encode_sub_20.rds')
meth[is.na(meth)] = 0
str(expr)
str(meth)

trainexpr <- expr[,1:11]
testexpr <- expr[,12:20]
trainmeth <- meth[,1:11]
testmeth <- meth[,12:20]

rm('expr')
rm('meth')


trainmodel <- function(trainexpr,testexpr,trainmeth) {
  library(data.table)
  library(fastcluster)
  
  trainmeth <- trainmeth[,colnames(trainexpr)]
  trainmeth[trainmeth==0] <- min(trainmeth[trainmeth>0])
  trainmeth[trainmeth==1] <- max(trainmeth[trainmeth<1])
  
  trainmeth <- log(trainmeth/(1-trainmeth))
  
  int <- intersect(rownames(trainexpr),rownames(testexpr))
  trainexpr <- trainexpr[int,]
  testexpr <- testexpr[int,]
  colnames(trainexpr) <- paste0('train_',1:ncol(trainexpr))
  colnames(testexpr) <- paste0('test_',1:ncol(testexpr))
  
  expr <- cbind(trainexpr,testexpr)
  
  qnem <- rowMeans(apply(expr,2,function(i) sort(i)))
  expr <- apply(expr,2,function(i) qnem[frank(i)])
  
  
  ### cross validation
  set.seed(12345)
  spid <- split(sample(colnames(expr)),ceiling(seq_along(colnames(expr))/(ncol(expr)/5)))
  perf <- NULL
  for (i in length(spid):1) {
    testid <- spid[[i]]
    trainid <- setdiff(colnames(expr),testid)
    cvexpr <- expr[,trainid]
    cvexpr <- cvexpr[rowSums(cvexpr) > 0,]
    cvexpr <- cvexpr[rowSums(cvexpr >= median(cvexpr[cvexpr>0])) >= 1,]
    
    m <- rowMeans(cvexpr)
    s <- sqrt((rowMeans(cvexpr*cvexpr) - m^2) * ncol(cvexpr)/(ncol(cvexpr)-1))
    cv <- s/m
    cvexpr <- cvexpr[cv >= quantile(cv,0.25),]
    
    qnem <- rowMeans(apply(cvexpr,2,function(i) sort(i)))
    gn <- row.names(cvexpr)
    cvexpr <- apply(cvexpr,2,function(i) qnem[frank(i)])
    row.names(cvexpr) <- gn
    cvtestexpr <- apply(expr[row.names(cvexpr),testid],2,function(i) qnem[frank(i)])
    row.names(cvtestexpr) <- gn
    testy <- meth[,testid]
    
    m <- rowMeans(cvexpr)
    s <- sqrt((rowMeans(cvexpr*cvexpr) - m^2) * ncol(cvexpr)/(ncol(cvexpr)-1))
    stde <- (cvexpr-m)/s
    
    for (cn in c(1000,2000,3000)) {
      print(cn)
      gclu <- kmeans(stde,cn)$cluster
      
      clucvexpr <- sapply(1:max(gclu),function(i) colMeans(cvexpr[names(gclu)[gclu==i],,drop=F]))
      clucvtestexpr <- sapply(1:max(gclu),function(i) colMeans(cvtestexpr[names(gclu)[gclu==i],,drop=F]))
      clucvexprmean <- colMeans(clucvexpr)
      clucvexprsd <- apply(clucvexpr,2,sd)
      
      trainx <- t((t(clucvexpr)-clucvexprmean)/clucvexprsd)
      testx <- t((t(clucvtestexpr)-clucvexprmean)/clucvexprsd)
      trainy <- meth[,trainid]
      
      trainymean <- rowMeans(trainy)
      trainysd <- apply(trainy,1,sd)
      trainy <- (trainy-trainymean)/trainysd
      
      crossx <- crossprod(trainx)
      
      for (lambda in 10^c(-2,-1,1,2)) {
        pred <- trainy %*% ((trainx %*% chol2inv(chol(crossx+ncol(trainy)*lambda*diag(ncol(trainx))))) %*% t(testx))
        pred <- pred * trainysd + trainymean
        perf <- rbind(perf,data.frame(cn=cn,lambda=lambda,mse=mean((pred[,1]-testy[,1])^2)))
      }
    }
  }
  print(perf) 
  perf <- aggregate(perf[,3],list(perf[,1],perf[,2]),mean)
  
  cn <- perf[which.min(perf[,3]),1]
  lambda <- perf[which.min(perf[,3]),2]
  
  expr <- expr[rowSums(expr) > 0,]
  expr <- expr[rowSums(expr >= median(expr[expr>0])) >= 1,]
  
  m <- rowMeans(expr)
  s <- sqrt((rowMeans(expr*expr) - m^2) * ncol(expr)/(ncol(expr)-1))
  cv <- s/m
  expr <- expr[cv >= quantile(cv,0.25),]
  qnem <- rowMeans(apply(expr,2,function(i) sort(i)))
  gn <- row.names(expr)
  expr <- apply(expr,2,function(i) qnem[frank(i)])
  row.names(expr) <- gn
  
  m <- rowMeans(expr)
  s <- sqrt((rowMeans(expr*expr) - m^2) * ncol(expr)/(ncol(expr)-1))
  stde <- (expr-m)/s
  
  gclu <- kmeans(stde,cn)$cluster
  
  cluexpr <- sapply(1:max(gclu),function(i) colMeans(expr[names(gclu)[gclu==i],,drop=F]))
  cluexprmean <- colMeans(cluexpr)
  cluexprsd <- apply(cluexpr,2,sd)
  trainx <- t((t(cluexpr)-cluexprmean)/cluexprsd)
  trainymean <- rowMeans(meth)
  trainysd <- apply(meth,1,sd)
  meth <- (meth-trainymean)/trainysd
  crossx <- crossprod(trainx)
  
  beta <- meth %*% (trainx %*% chol2inv(chol(crossx+ncol(trainy)*lambda*diag(ncol(trainx)))))
  names(gclu) <- sub('\\..*','',names(gclu))
  row.names(expr) <- sub('\\..*','',row.names(expr))
  list(beta=beta,trainexpr=expr,genecluster=gclu,trainxmean=cluexprmean,trainxsd=cluexprsd,trainymean=unname(trainymean),trainysd=unname(trainysd))
}



