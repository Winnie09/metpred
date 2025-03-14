
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
  # trainmeth[trainmeth==0] <- min(trainmeth[trainmeth>0]) 
  # trainmeth[trainmeth==1] <- max(trainmeth[trainmeth<1]) 
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

