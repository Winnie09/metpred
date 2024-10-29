## Task: Emily's mouse pancreas data, cross validation.
## Prediction step. 
## train: celltype-specific gene expression (gene by sample:celltype, celltypes are from NNLS)
## and cell-type-primary DNAm (CpG by sample:celltype, celltypes are from LCM)
## test:celltype-specific gene expression
## result: predicted celltype-primary DNAm (CpG by sample:celltype)
i <- as.character(commandArgs(trailingOnly = T))
print(i)

setwd('/home/whou10/data/whou/metpred/')
rdir <- paste0('evaluate/mousepancreas/celltypespecific_cv/res/',i)
dir.create(rdir, showWarnings = F, recursive = T)
# i = 'A46_F_P_K4_D4'
r <- readRDS('data/mousepancreas/spatial/proc/all/cse.rds')
w <- readRDS('data/mousepancreas/wgbs/bs.rds')
# > str(r)
#  num [1:19465, 1:280] 0.000759 0 0 0.043464 0.007568 ...
#  - attr(*, "dimnames")=List of 2
#   ..$ : chr [1:19465] "Xkr4" "Rp1" "Sox17" "Lypla1" ...
#   ..$ : chr [1:280] "A29_M_P_K4_D2:acinar" "A29_M_P_K4_D2:adipose" "A29_M_P_K4_D2:adm" "A29_M_P_K4_D2:b_cell" ...
# > str(w)
#  num [1:19634603, 1:41] 0.767 0.764 0.764 0.759 0.753 ...
#  - attr(*, "dimnames")=List of 2
#   ..$ : chr [1:19634603] "chr1:3000827" "chr1:3001007" "chr1:3001018" "chr1:3001277" ...
#   ..$ : chr [1:41] "A29_M_P_K4_D2:adm" "A29_M_P_K4_D2:duct" "A31_M_P_K4_D4:adm" "A31_M_P_K4_D4:duct" ...
int <- intersect(colnames(r),colnames(w))
samp <- sub(':.*','',int)
# > str(samp)
#  chr [1:41] "A29_M_P_K4_D2" "A29_M_P_K4_D2" "A34_M_P_K4_D2" "A34_M_P_K4_D2" ...
trainid <- int[which(samp != i)]
testid <- setdiff(int, trainid)
source('software/trainpredict.R')
source('/home/whou10/scratch16/whou10/resource/startup.R')

## group the CpGs from its standard deviation in testing samples smallest to largets slots
## run the prediction for each CpG group
rs <- rowsds(w[,testid])
rs2 <- cut(rs,seq(0,1,0.05))
names(rs2) <- names(rs)

for (rs2.s in unique(rs2)){
  fn <- paste0(rdir,'/testsd_cpg_', which(levels(rs2) == rs2.s), '.rds')
  print(rs2.s)
  if (file.exists(fn)){
    next
  } 
  pred <- trainpredict(trainexpr=r[,trainid],testexpr=r[,testid],trainmeth=w[names(rs2)[which(rs2 == rs2.s)], trainid])
  saveRDS(pred,file=fn)
}

