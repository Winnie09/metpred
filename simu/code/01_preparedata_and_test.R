nSample = as.numeric(commandArgs(trailingOnly = T)[1]) ## number of sample in each group, 2, 4
nCell = as.numeric(commandArgs(trailingOnly = T)[2]) ## number of cells in each sample:2 4 10 100 200
propCpG = as.numeric(commandArgs(trailingOnly = T)[3]) ## 0.2, 0.3, 0.4
probRead= as.numeric(commandArgs(trailingOnly = T)[4]) ## 0.3 0.4 0.5
# nSample = 2
# nCell = 2
# propCpG = 0.2
# probRead= 0.3
rdir <- paste0('/home-4/whou10@jhu.edu/work-zfs/whou10/metpred/simu/res/', nSample, '_', nCell, '_', propCpG, '_', probRead, '/')
dir.create(rdir, recursive = TRUE)
datafile <- paste0('/home-4/whou10@jhu.edu/work-zfs/whou10/metpred/simu/data/', nSample, '_', nCell, '_', propCpG, '_', probRead,'.rds')

logistic <- function(x) exp(x)/(1+exp(x))
logit <- function(x) log(x/(1-x))
source('/home-4/whou10@jhu.edu/work-zfs/whou10/trajectory_variability/function/01_function.R')
source('/home-4/whou10@jhu.edu/scratch/Wenpin/raisin/software/raisin/raisin.R')
source('/home-4/whou10@jhu.edu/work-zfs/whou10/metpred/software/findDMR.R')
func <- func[grepl('saver',names(func))]
func <- func[!names(func) %in% 'limmacell_saver']

fullmat <- readRDS(paste0('/home-4/whou10@jhu.edu/work-zfs/whou10/metpred/nullsimu/data/data/', nSample,'To', nSample, '_allcell/pred_ery.rds'))
ap <- as.numeric(sub('BM','',sub(':.*', '', colnames(fullmat))))
as <- sub('_.*','', colnames(fullmat))
ind <- NULL

Res <- NULL

mat <- lapply(unique(as), function(i){
  # set.seed(9884)
  set.seed(1)
  #set.seed(6944)
  id <- sample(which(as == i), nCell, replace = TRUE)
  fullmat[, id]
})
mat <- do.call(cbind, mat)
colnames(mat) <- paste0(colnames(mat), '_', 1:ncol(mat))

## add signals
set.seed(12345)
addgene <- sample(rownames(mat),propCpG*nrow(mat))
cellid <- which(as.numeric(sub('BM','',sub(':.*','',colnames(mat)))) %in% c(1,2,5,6))
addmat <- mat[sample(setdiff(rownames(mat),rownames(addgene)),length(addgene)),sample(colnames(mat),length(cellid))]
mat[addgene, cellid] <- logistic(logit(mat[addgene, cellid]) + logit(addmat * probRead))
rownames(mat)[rownames(mat) %in% addgene] <- paste0('add:',rownames(mat)[rownames(mat) %in% addgene])

## prepare for test
sample <- sub(':.*','',colnames(mat))
samplename <- unique(sample)
sum(as.numeric(sub('BM','',samplename)) %in% c(1,2,5,6) + 0 == 1) < length(samplename)
group <- as.numeric(sub('BM','',samplename)) %in% c(1,2,5,6) + 1
expr <- logit(mat) ## transform to -inf, +inf

## test our method
design <- data.frame(sample=samplename,feature=group)
res1 <- RAISINtest(RAISINfit(expr,sample,design,testtype='unpaired',filtergene=F)) ## -inf, +inf
saveRDS(res1, file=paste0(rdir, 'ourmethod.rds'))

## test other methods
ind <- NULL
expr <- mat
for (metid in 1:length(func)) {
  fun <- func[[metid]]
  res <- fun()
  saveRDS(res,file=paste0(rdir, names(func)[metid],'.rds'))  
}



