## debug: why acrosssample pcc very low??
library(ggplot2)
pal <- c('Ramp'='red','Nearest gene'='orange','Permute'='royalblue')
pd <- NULL
type <- 'sample'
print(type)
rdir <- '/home/whou10/data/whou/metpred/evaluate/pbmc/predsc/imp/perf/'
p1 <- data.frame(method='Nearest gene',type=type,readRDS(paste0(rdir, 'neargene/pred_var_10000_CpG_across',type,'_pcc.rds')),stringsAsFactors = F)
p2 <- data.frame(method='Permute',type=type,readRDS(paste0(rdir,'permu/pred_var_10000_CpG_across',type,'_pcc.rds')),stringsAsFactors = F)
p3 <- data.frame(method='Ramp',type=type,readRDS(paste0(rdir, 'ramp/pred_var_10000_CpG_across',type,'_pcc.rds')),stringsAsFactors = F)
pd <- rbind(p1,p2,p3)
pd$sd <- as.factor(as.character(pd$sd))
str(pd)


tapply(p3$cor, p3$sd, median, na.rm=T)
tmp = p3[p3$sd=='(0.35,0.4]', ]
tmp = tmp[tmp$cor < median(tmp$cor, na.rm=T), ]
tmp = tmp[order(tmp$cor, decreasing = T),]


## load gs and pred
pred <- readRDS('/home/whou10/data/whou/metpred/evaluate/pbmc/predsc/imp/res/eval/eval_pred_ctmean.rds')
me <- readRDS('/home/whou10/data/whou/metpred/evaluate/pbmc/predsc/imp/res/eval/eval_goldstandard_ctmean.rds')

## check cpg
for (i in 1:20){
  print(x <- pred[rownames(tmp)[i], ])
  print(y <- me[rownames(tmp)[i], ])
  print(cor(x, y))  
}
selcpg = rownames(tmp)[4] ## 'chr16:74407150'
# cd56_nk cd14_monocytes        b_cells        t_cells 
# 0.8619493      0.2922073      0.5192954      0.6282576 
# cd56_nk cd14_monocytes        b_cells        t_cells 
# 0.01192019     0.73841685     0.86443540     0.50343602 
# [1] -0.8519115

## load training methylation
w.tr = readRDS('/home/whou10/data/whou/metpred/data/pbmc/imp/res/wgbs_train.rds')
str(w.tr)
# # num [1:30673446, 1:66] 0.744 0.743 0.741 0.741 0.74 ...
# # - attr(*, "dimnames")=List of 2
# # ..$ : chr [1:30673446] "chr1_10469" "chr1_10471" "chr1_10484" "chr1_10489" ...
# # ..$ : chr [1:66] "blueprint_venous_blood-C000S5-CD14-positive_CD16-negative_classical_monocyte" "blueprint_venous_blood-C000S5-mature_neutrophil" "blueprint_venous_blood-C0010K-CD14-positive_CD16-negative_classical_monocyte" "blueprint_venous_blood-C0010K-mature_neutrophil" ...
# rs2 = readRDS(paste0(ddir, 'wgbs_train_nonNA_CpG_rowsds_cut.rds'))
# w.tr = w.tr[names(rs2)[rs2 == levels(rs2)[cpggroup]], , drop = F]
# if (nrow(w.tr) < 1){
#   stop('No CpG in this rowsds group!')
# }
w.tr[selcpg, grep('CD56', colnames(w.tr))]
x = w.tr[,"blueprint_cord_blood-C0067N-cytotoxic_CD56-dim_natural_killer_cell"]
y = me[,'cd56_nk']
int = intersect(names(x), names(y))
x = x[int]
y = y[int]
cor(x,y)
v = x-y
summary(v)


## identify CpGs that have large difference between training and testing cell types
w.te = readRDS('/home/whou10/data/whou/metpred/data/pbmc/imp/res/wgbs_test.rds')
# > colnames(w.te)
# [1] "blueprint_cord_blood-S01E8O-cytotoxic_CD56-dim_natural_killer_cell"          
# [2] "blueprint_venous_blood-C004SQ-CD14-positive_CD16-negative_classical_monocyte"
# [3] "blueprint_cord_blood-C005PS-CD14-positive_CD16-negative_classical_monocyte"  
# [4] "blueprint_cord_blood-S000RD-CD14-positive_CD16-negative_classical_monocyte"  
# [5] "blueprint_tonsil-T14_11-germinal_center_B_cell"                              
# [6] "blueprint_venous_blood-S009W4-CD4-positive_alpha-beta_T_cell"                
# [7] "blueprint_venous_blood-S014QS-central_memory_CD4-positive_alpha-beta_T_cell" 
# [8] "blueprint_venous_blood-S014QS-effector_memory_CD4-positive_alpha-beta_T_cell"
# [9] "blueprint_cord_blood-S00C2F-CD8-positive_alpha-beta_T_cell"   

w.te2 = w.te[, -c(1,3,4,9)]
## testing cell type: similar to pbmc cell type; similar to training celltype methylation
## testing cell type:
# c('b_cells', 'cd14_monocytes', 'cd56_nk', 'cytotoxic_t', 'memory_t', 'naive_t','regulatory_t')

pred.agg = readRDS('/home/whou10/data/whou/metpred/evaluate/pbmc/predsc/imp/res/pb/pred_allcpg.rds')
w.te3 = w.te2[rownames(pred.agg),]
w.te.agg = aggregatefunc2(d = w.te3, 
                          by = c('cd14_monocytes', 'b_cells', rep('t_cells', 3)), 
                          fun = 'mean')
pred.agg = pred.agg[, colnames(w.te.agg)]

pd = acrossRowCor_plotdata(pred = pred.agg, goldstandard = w.te.agg)
saveRDS(pd, 'tmp_ramp_pd.rds')
