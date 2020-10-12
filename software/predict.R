predict <- function(expr,model,impute=F) {
  # impute: impute the genes that are in training set but not in prediction set.
  library(data.table)
  library(SMUT)

  # Quantile normalization testing data to training
  row.names(expr) <- sub('\\..*','',row.names(expr))  
  intgene <- intersect(row.names(expr),row.names(model$trainexpr))
  qnem <- rowMeans(apply(model$trainexpr[intgene,],2,function(i) sort(i)))
  expr <- expr[intgene,]
  gn <- row.names(expr)
  expr <- apply(expr,2,function(i) qnem[frank(i,ties.method='min')])
  row.names(expr) <- gn
  
  # Do we want to impute genes that are in the training set but not in testing set?
  if (!impute) {
    cluexpr <- sapply(1:max(model$genecluster),function(i) colMeans(expr[intersect(row.names(expr),names(model$genecluster)[model$genecluster==i]),,drop=F]))
  } else {
    trainx <- model$trainexpr[row.names(expr),]
    if (ncol(trainx) > 1000) trainx <- trainx[,sample(1:ncol(trainx),1000)]
    m <- rowMeans(trainx)
    s <- sqrt((rowMeans(trainx*trainx) - m^2) * ncol(trainx)/(ncol(trainx)-1))
    stde <- (trainx-m)/s
    print('begin kmeans')
    trainxclu <- kmeans(stde,nrow(stde)/10,iter.max = 10000)$cluster
    trainxcluexpr <- sapply(1:max(trainxclu),function(i) colMeans(trainx[names(trainxclu)[trainxclu==i],,drop=F]))
    trainy <- model$trainexpr[setdiff(row.names(model$trainexpr),row.names(expr)),]
    
    print('begin beta')
    lambda <- 1
    logtrainy <- trainy
    if (min(logtrainy)==0) logtrainy[logtrainy==0] <- min(logtrainy[logtrainy>0])
    logtrainy <- log(logtrainy)
    beta <- eigenMapMatMult(logtrainy , eigenMapMatMult(trainxcluexpr , chol2inv(chol(crossprod(trainxcluexpr) + lambda * nrow(trainxcluexpr) * diag(ncol(trainxcluexpr))))))
    
    print('begin pred')
    predxcluexpr <- sapply(1:max(trainxclu),function(i) colMeans(expr[names(trainxclu)[trainxclu==i],,drop=F]))
    predy <- exp(eigenMapMatMult(beta, t(predxcluexpr)))
    
    finalexpr <- matrix(0,nrow=nrow(model$trainexpr),ncol=ncol(expr),dimnames=list(row.names(model$trainexpr),colnames(expr)))
    finalexpr[row.names(expr),colnames(expr)] <- expr
    finalexpr[row.names(predy),colnames(predy)] <- predy
    
    qnem <- rowMeans(apply(model$trainexpr,2,function(i) sort(i)))
    dn <- dimnames(finalexpr)
    finalexpr <- apply(finalexpr,2,function(i) qnem[frank(i)])
    dimnames(finalexpr) <- dn
    
    cluexpr <- sapply(1:max(model$genecluster),function(i) colMeans(finalexpr[names(model$genecluster)[model$genecluster==i],,drop=F]))
  }
  
  # Predict the standardized value and add back mean and sd
  testx <- t((t(cluexpr)-model$trainxmean)/model$trainxsd)
  testx[is.na(testx)] <- 0
  pred <- eigenMapMatMult(model$beta, t(testx))
  pred <- pred * model$trainysd + model$trainymean
  pred <- exp(pred)
  pred <- pred/(1+pred)
  row.names(pred) <- row.names(model$beta)
  colnames(pred) <- colnames(expr)
  pred
}

