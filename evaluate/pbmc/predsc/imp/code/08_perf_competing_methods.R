## calculate pbmc perf
rm(list=ls())
method = commandArgs(trailingOnly = T)[1] ##. neargene or permu
key = 'pred_var_10000_CpG'
print(method)

# key = 'pred_varmedian_0673446_CpG'

source('/home/whou10/scratch16/whou10/resource/myfunc/01_function.R')
ddir <- '/home/whou10/data/whou/metpred/data/pbmc/imp/res/'
w.te = readRDS(paste0(ddir, 'wgbs_test.rds'))
# > str(w.te)
# num [1:30673446, 1:9] 0.772 0.771 0.767 0.766 0.765 ...
# - attr(*, "dimnames")=List of 2
# ..$ : chr [1:30673446] "chr1_10469" "chr1_10471" "chr1_10484" "chr1_10489" ...
# ..$ : chr [1:9] "blueprint_cord_blood-S01E8O-cytotoxic_CD56-dim_natural_killer_cell" "blueprint_venous_blood-C004SQ-CD14-positive_CD16-negative_classical_monocyte" "blueprint_cord_blood-C005PS-CD14-positive_CD16-negative_classical_monocyte" "blueprint_cord_blood-S000RD-CD14-positive_CD16-negative_classical_monocyte" ...

rdir = '/home/whou10/data/whou/metpred/evaluate/pbmc/predsc/imp/'

## load prediction
res = readRDS(paste0(rdir, '/competing/combine/', method, '.rds'))
# > str(res)
# num [1:30512213, 1:4] 0.712 0.712 0.712 0.712 0.712 ...
# - attr(*, "dimnames")=List of 2
# ..$ : Named chr [1:30512213] "chr1_10469" "chr1_10471" "chr1_10484" "chr1_10489" ...
# .. ..- attr(*, "names")= chr [1:30512213] "chr1:10469" "chr1:10471" "chr1:10484" "chr1:10489" ...
# ..$ : chr [1:4] "b_cells" "cd14_monocytes" "cd56_nk" "t_cells"


## load filtered-in CpGs
ramp = readRDS(paste0(rdir, 'res/pred/pred_var_100000_CpG.rds'))
# > str(pred)
# num [1:100000, 1:37740] 6.49e-11 2.07e-67 1.82e-65 1.98e-61 4.04e-67 ...
# - attr(*, "dimnames")=List of 2
# ..$ : chr [1:100000] "chr14_34541540" "chr16_22435056" "chr16_22435145" "chr16_22435227" ...
# ..$ : Named chr [1:37740] "b_cells:AAACATACAATGCC-1" "b_cells:AAACATACACGCAT-1" "b_cells:AAACATACGAATAG-1" "b_cells:AAACATACGTGTCA-1" ...
# .. ..- attr(*, "names")= chr [1:37740] "test_1" "test_2" "test_3" "test_4" 

## prepare predition and goldstandard: filter-in CpG by bulk cell type
o = intersect(rownames(ramp), rownames(res))
pred.agg = res[o, ]
w.te2 = w.te[o,]
w.te.agg = aggregatefunc2(d = w.te2, 
                          by = c('cd56_nk', rep('cd14_monocytes', 3), 'b_cells', rep('t_cells', 4)), 
                          fun = 'mean')
pred.agg = pred.agg[, colnames(w.te.agg)]

## across-sample PCC
pd = acrossRowCor_plotdata(pred = pred.agg, goldstandard = w.te.agg)
saveRDS(pd, paste0(rdir, 'perf/', method, '/', key, '_acrosssample_pcc.rds'))


## across-cpg PCC
id = complete.cases(pred.agg) & complete.cases(w.te.agg)
pd2 = acrossRowCor_plotdata(pred = t(pred.agg[id,]), 
                            goldstandard = t(w.te.agg[id,]))
saveRDS(pd2, paste0(rdir, 'perf/', method, '/', key, '_acrosscpg_pcc.rds'))

