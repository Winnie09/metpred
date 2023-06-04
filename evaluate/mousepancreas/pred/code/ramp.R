i <- as.character(commandArgs(trailingOnly = T))

if (file.exists(paste0('/home/whou10/data/whou/metpred/evaluate/mousepancreas/pred/ramp/',i,'.rds'))){
  stop('File exists. ')
}

dir.create(paste0('/home/whou10/data/whou/metpred/evaluate/mousepancreas/pred/ramp/',i))


r <- readRDS('/home/whou10/data/whou/metpred/data/mousepancreas/spatial/cse.rds')
w <- readRDS('/home/whou10/data/whou/metpred/data/mousepancreas/wgbs/bs.rds')
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
library(Matrix)
e <- readRDS('/home/whou10/data/whou/metpred/data/mousepancreas/spatial/normexpr.rds')
se <- sub(':.*','',colnames(e))
e <- e[,se == i]

# > table(se)
# se
#  A29_M_P_K4_D2  A31_M_P_K4_D4  A34_M_P_K4_D2  A46_F_P_K4_D4  A47_F_P_K4_D7 
#           2764           3102           2385           4124           3079 
#    A53_F_LM_D2 A60_M_WT_LM_D2  A66_M_P_K4_D7    A72_F_LM_D7    A73_F_LM_D4 
#           2979           2783           3112           3270           3742 
#          D2FC1          D2FP1          D2MC1          D2MP1          D4FC3 
#           1849           1400           1030            788           1697 
#          D4FP4          D4MC4          D7FC6          D7FP6          D7MC6 
#           1455           1563           2105           2053           2139 
          
source('/home/whou10/data/whou/metpred/software/trainpredict.R')
source('/home/whou10/scratch16/whou10/resource/startup.R')

## group the CpGs from its standard deviation in testing samples smallest to largets slots
## run the prediction for each CpG group
rs <- rowsds(w[,testid])
rs2 <- cut(rs,seq(0,1,0.05))
names(rs2) <- names(rs)

for (rs2.s in unique(rs2)){
  print(rs2.s)
  pred <- trainpredict(trainexpr=r[,trainid],testexpr=e,trainmeth=w[names(rs2)[which(rs2 == rs2.s)], trainid])
  saveRDS(pred,file=paste0('/home/whou10/data/whou/metpred/evaluate/mousepancreas/pred/ramp/',i,'/testsd_cpg_', which(levels(rs2) == rs2.s), '.rds'))
}


