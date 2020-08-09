library(SMUT)
trainmodel <- function(expr,meth) {
  set.seed(12345)
  library(data.table)
  library(fastcluster)
  
  # avoid 0 and 1 to enter logit function
  meth[meth==0] <- min(meth[meth>0])
  meth[meth==1] <- max(meth[meth<1])
  meth <- log(meth/(1-meth))
  
  # Filtering low expressed genes and quantile normalization
  expr <- expr[rowMeans(expr > 0) >= 0.1,]
  qnem <- rowMeans(apply(expr,2,function(i) sort(i)))
  gn <- row.names(expr)
  expr <- apply(expr,2,function(i) qnem[frank(i,ties.method='min')])
  row.names(expr) <- gn
  
  # Filtering low variable genes
  m <- rowMeans(expr)
  s <- sqrt((rowMeans(expr*expr) - m^2) * ncol(expr)/(ncol(expr)-1))
  mod <- loess(s~m)
  expr <- expr[resid(mod) > 0,]
  
  # Gene standardization and clustering
  m <- rowMeans(expr)
  s <- sqrt((rowMeans(expr*expr) - m^2) * ncol(expr)/(ncol(expr)-1))
  stde <- (expr-m)/s
  gclu <- kmeans(stde,nrow(stde)/10,iter.max=10000)$cluster
  print('kmeans done')
  
  ### cross validation to select lambda
  spid <- split(sample(colnames(expr)),ceiling(seq_along(colnames(expr))/(ncol(expr)/5)))
  llist <- c(0.1,1,10,100,1000,10000)
  perf <- matrix(0,nrow=nrow(meth),ncol=length(llist),dimnames=list(row.names(meth),llist))
  for (i in length(spid):1) {
    print(i)
    testid <- spid[[i]]
    trainid <- setdiff(colnames(expr),testid)
    cvexpr <- expr[,trainid]
    cvtestexpr <- expr[,testid]
    testy <- meth[,testid]
    
    # Get averaged gene expression in each cluster for training and cv
    clucvexpr <- sapply(1:max(gclu),function(i) colMeans(cvexpr[names(gclu)[gclu==i],,drop=F]))
    clucvtestexpr <- sapply(1:max(gclu),function(i) colMeans(cvtestexpr[names(gclu)[gclu==i],,drop=F]))
    clucvexprmean <- colMeans(clucvexpr)
    clucvexprsd <- apply(clucvexpr,2,sd)
    
    # Standardize using mean and sd in training set
    trainx <- t((t(clucvexpr)-clucvexprmean)/clucvexprsd)
    testx <- t((t(clucvtestexpr)-clucvexprmean)/clucvexprsd)
    trainy <- meth[,trainid]
    
    trainymean <- rowMeans(trainy)
    trainysd <- apply(trainy,1,sd)
    trainy <- (trainy-trainymean)/trainysd
    
    crossx <- crossprod(trainx)
    
    #    if (nrow(trainx) > ncol(trainx)) llist <- c(0,llist)
    # Fit ridge regression, calculate MSE
    for (lambda in llist) {
      pred <- eigenMapMatMult(trainy , eigenMapMatMult(eigenMapMatMult(trainx, chol2inv(chol(crossx+ncol(trainy)*lambda*diag(ncol(trainx))))) , t(testx)))
      pred <- pred * trainysd + trainymean
      mse <- rowMeans((pred-testy)^2)
      mse[is.na(mse)] <- 0
      perf[,as.character(lambda)] <- perf[,as.character(lambda)]+mse
    }
  }
  rm('trainy')
  # Choose best lambda
  lambda <- as.numeric(colnames(perf)[apply(perf,1,which.min)])
  
  # With the best lambda, retrain using all data (no cv)
  cluexpr <- sapply(1:max(gclu),function(i) colMeans(expr[names(gclu)[gclu==i],,drop=F]))
  cluexprmean <- colMeans(cluexpr)
  cluexprsd <- apply(cluexpr,2,sd)
  trainx <- t((t(cluexpr)-cluexprmean)/cluexprsd)
  trainymean <- rowMeans(meth)
  trainysd <- apply(meth,1,sd)
  meth <- (meth-trainymean)/trainysd
  crossx <- crossprod(trainx)
  
  beta <- matrix(0,nrow=nrow(meth),ncol=ncol(trainx),dimnames=list(row.names(meth),colnames(trainx)))
  for (l in unique(lambda)) {
    id <- which(lambda==l)
    beta[id,] <- eigenMapMatMult(meth[id,,drop=F], eigenMapMatMult(trainx , chol2inv(chol(crossx+ncol(meth)*l*diag(ncol(trainx))))))  
  }
  names(gclu) <- sub('\\..*','',names(gclu))
  row.names(expr) <- sub('\\..*','',row.names(expr))
  list(beta=beta,trainexpr=expr,genecluster=gclu,trainxmean=cluexprmean,trainxsd=cluexprsd,trainymean=unname(trainymean),trainysd=unname(trainysd),lambda=lambda)
}
