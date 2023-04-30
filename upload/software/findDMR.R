## identify DMR
func <- list()
library(parallel)

func[['limmacell_saver']] <- function(expr.=expr,sample.=sample,samplename.=samplename,group.=group,ind.=ind) {
  suppressMessages(library(limma))
  if (is.null(ind)) {
    design <- cbind(Grp1=1,Grp2vs1=as.numeric(sample %in% samplename[group==1]))  
  } else {
    design <- cbind(Grp1=1,Grp2vs1=as.numeric(sample %in% samplename[group==1]),ind=model.matrix(~x-1,data.frame(x=ind[sample])))
  }
  options(digits=3)
  block <- as.numeric(as.factor(sample))
  dupcor <- duplicateCorrelation(expr,design,block=block)
  if (is.null(ind)) {
        fit <- lmFit(expr, design, block = block, correlation = dupcor$consensus)
  } else {
  tryCatch(fit <- lmFit(expr, design, block = block, correlation = dupcor$consensus), error=function(e) {})
  if (!exists('fit')) {
        design <- cbind(Grp1=1,Grp2vs1=as.numeric(sample %in% samplename[group==1]))
        block <- as.numeric(as.factor(sample))
        dupcor <- duplicateCorrelation(expr,design,block=block)
        fit <- lmFit(expr, design, block = block, correlation = dupcor$consensus)
  }
  }
  fit <- eBayes(fit)
  res <- topTable(fit,coef=2,number=nrow(expr))
  colnames(res)[colnames(res)=='P.Value'] <- "pvalue"
  colnames(res)[colnames(res)=='t'] <- "stat"
  colnames(res)[which(colnames(res) == 'adj.P.Val')] <- 'FDR'
  res[order(res[,'pvalue'],-abs(res[,'stat'])),]
}

func[['limma_saver']] <- function(expr.=expr,sample.=sample,samplename.=samplename,group.=group,ind.=ind) {
  suppressMessages(library(limma))
  if (is.null(ind) || sum(duplicated(ind))==0) {
    design <- cbind(Grp1=1,Grp2vs1=as.numeric(group==group[1]))
  } else {
    design <- cbind(Grp1=1,Grp2vs1=as.numeric(group==group[1]),ind=model.matrix(~x-1,data.frame(x=ind)))
  }
  sampsum <- sapply(samplename,function(us) rowMeans(expr[,sample==us,drop=F]))
  options(digits=3)
  fit <- lmFit(sampsum, design=design)
  fit <- eBayes(fit)
  res <- topTable(fit,coef=2,number=nrow(expr))
  colnames(res)[colnames(res)=='P.Value'] <- "pvalue"
  colnames(res)[colnames(res)=='t'] <- "stat"
  colnames(res)[which(colnames(res) == 'adj.P.Val')] <- 'FDR'
  res[order(res[,'pvalue'],-abs(res[,'stat'])),]
}

func[['MAST_saver']] <- function(expr.=expr,sample.=sample,samplename.=samplename,group.=group) {
  suppressMessages(library(MAST))
  library(reshape2)
  # did not scale cdr <- scale(colMeans(count > 0))
  cluster <- rep('clu1',ncol(expr))
  cluster[sample %in% samplename[group==2]] <- 'clu2'
  sca <- FromMatrix(expr,data.frame(wellKey=colnames(expr),cluster=cluster), data.frame(primerid=row.names(expr),Gene=row.names(expr)))
  zlmCond <- zlm(~cluster, sca)
  summaryCond <- summary(zlmCond, doLRT="clusterclu2")
  summaryDt <- as.data.frame(summaryCond$datatable)
  pval <- summaryDt[summaryDt$component == "H" & summaryDt$contrast == "clusterclu2",c(1,4)]
  lfc <- summaryDt[summaryDt$component == "logFC" & summaryDt$contrast == "clusterclu2",c(1,7)]
  combine <- merge(pval,lfc)
  colnames(combine) <- c("Gene","pvalue","Log-foldchange")
  row.names(combine) <- combine[,1]
  res <- combine
  colnames(res)[colnames(res)=='Log-foldchange'] <- "stat"
  res$FDR <- p.adjust(res$pvalue, method = 'fdr')
  res[order(res[,'pvalue'],-abs(res[,'stat'])),]
}


func[['scDD_saver']] <- function(expr.=expr, sample.=sample,samplename.=samplename,group.=group) {
  suppressMessages(library(scDD))
  suppressMessages(library(SingleCellExperiment))
  condition <- as.numeric(sample %in% samplename[group==2]) + 1
  names(condition) <- colnames(expr)
  d <- SingleCellExperiment(assays=list(normcounts=expr),colData=data.frame(condition))
  prior_param=list(alpha=0.01, mu0=0, s0=0.01, a0=0.01, b0=0.01)
  res <- scDD(d, prior_param=prior_param, testZeroes=T,param=BiocParallel::MulticoreParam(workers=2),categorize=FALSE)
  res <- res@metadata$Genes
  colnames(res)[colnames(res)=='combined.pvalue'] <- 'pvalue'
  res$FDR <- p.adjust(res$pvalue, method = 'fdr')
  res[order(res[,'pvalue']),]
  
}

func[['t_saver']] <- function(expr.=expr, sample.=sample,samplename.=samplename,group.=group) {
  id1 <- which(sample %in% samplename[group==1])
  id2 <- which(sample %in% samplename[group==2])
  pval <- t(apply(expr,1,function(i) {
    if (length(unique(i[id1])) == 1 & length(unique(i[id2])) == 1) {
      c(1,0)
    } else {
      tmp <- t.test(i[id1],i[id2])
      c(tmp$p.value,tmp$statistic)
    }
  }))
  fdr <- p.adjust(pval[,1],method='fdr')
  res <- data.frame(Gene=row.names(pval),pvalue=pval[,1],FDR=fdr,stat=pval[,2],stringsAsFactors = F)
  res[order(res[,'pvalue'],-abs(res[,'stat'])),]
}

func[['wilcoxon_saver']] <- function(expr.=expr,sample.=sample,samplename.=samplename,group.=group) {
  id1 <- which(sample %in% samplename[group==1])
  id2 <- which(sample %in% samplename[group==2])
  pval <- t(apply(expr,1,function(i) {
    tmp <- wilcox.test(i[id1],i[id2])
    c(tmp$p.value,tmp$statistic)
  }))
  fdr <- p.adjust(pval[,1],method='fdr')
  res <- data.frame(Gene=row.names(pval),pvalue=pval[,1],FDR=fdr,stat=pval[,2],stringsAsFactors = F)
  res[order(res[,'pvalue'],-abs(res[,'stat'])),]
}


func[['glmm_saver']] <- function(expr.=expr,sample.=sample,samplename.=samplename,group.=group,ind.=ind) {
library(lmerTest)
 res <- mclapply(1:nrow(expr),function(i) {
    if (i %% 100 == 0) print(i)
    if (is.null(ind)) {
      tryCatch(suppressMessages(suppressWarnings(fit <- lmer(expr ~ group + (1|sample),data=data.frame(expr=expr[i,],group=group[match(sample,samplename)],sample=sample)))),error=function(e) {})
    } else {
      tryCatch(suppressMessages(suppressWarnings(fit <- lmer(expr ~ group + (group|ind),data=data.frame(expr=expr[i,],group=group[match(sample,samplename)],sample=sample,ind=ind[sample])))),error=function(e) {})
    }
    if (exists('fit')) {
    summary(fit)$coefficients[2,c(4,5)]
    } else {
      c(0,1)
    }
  },mc.cores=detectCores()-1)
  res <- do.call(rbind,res)
  res <- data.frame(res)
  colnames(res) <- c('stat','pvalue')
  row.names(res) <- row.names(expr)
  res$FDR <- p.adjust(res[,'pvalue'],method='fdr')
  res[order(res[,'pvalue'],-abs(res[,'stat'])),]
}

func[['BSmooth_saver']] <- function(expr.=expr, sample.=sample,samplename.=samplename,group.=group) {
  sampsum <- sapply(samplename,function(us) rowMeans(expr[,sample==us,drop=F]))
  id1 <- which(group==1)
  id2 <- which(group==2)
  pval <- t(apply(sampsum,1,function(i) {
    if (length(unique(i[id1])) == 1 & length(unique(i[id2])) == 1) {
      c(1,0)
    } else {
      tmp <- t.test(i[id1],i[id2])
      c(tmp$p.value,tmp$statistic)
    }
  }))
  fdr <- p.adjust(pval[,1],method='fdr')
  res <- data.frame(Gene=row.names(pval),pvalue=pval[,1],FDR=fdr,stat=pval[,2],stringsAsFactors = F)
  res[order(res[,'pvalue'],-abs(res[,'stat'])),]
}

