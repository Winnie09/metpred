ddir <- '/home/whou10/scratch16/whou10/trajectory_variability/hca_bone_marrow_data_analysis/real/build_from_tree_variability/result/erythroid/'
expr = readRDS(paste0(ddir, 'input_expr.rds'))
rownames(expr) = sub('.*:','',rownames(expr))

ge <- readRDS('/home/whou10/data/whou/metpred/data/encode_wgbs/processed/GRCh38_normal/rna.rds')
me <- readRDS('/home/whou10/data/whou/metpred/data/encode_wgbs/processed/GRCh38_normal/wgbs.rds')

# new program
source('/home/whou10/data/whou/metpred/software/trainpredict.R')
pred <- trainpredict(trainexpr=ge,testexpr=expr,trainmeth=me[1:1e5,])   
diff = readRDS('/home/whou10/data/whou/metpred/evaluate/hsc/erythroid/selcpg/mep_hsc/meandiff.rds')
pred2 = pred[names(sort(diff[rownames(pred)], decreasing = T)), ]
saveRDS(pred2,file='/home/whou10/data/whou/metpred/evaluate/hsc/erythroid/res/pred.rds')





library(data.table)
d <- fread('/home/whou10/data/whou/metpred/data/hsc/raw/rna/GSE87195_rnaseq_ensT_all.csv.gz',data.table=F)
rownames(d) <- d[,1]
d <- as.matrix(d[,-1])
w <- readRDS('/home/whou10/data/whou/metpred/data/hsc/proc/wgbs.rds')

rnact <- paste0(sub('.*_','',sub('_100','',colnames(d))),'_',sub('_.*','',sub('RNA_','',colnames(d))))
wct <- sub('_R[12].*','',sub('_1000','',sub('GSM[0-9]*_','',colnames(w))))
colnames(d) <- rnact
colnames(w) <- wct

int <- intersect(rnact,wct)
d <- d[,int]
w <- w[,int]

gtf <- fread('/home/whou10/scratch16/whou10/resource/gtf/grch38.gtf',data.table=F)
gtf <- gtf[gtf[,3]=='transcript',]
trid <- sub('\\..*','',sub('.*transcript_id "','',gtf[,9]))
gn <- sub('".*','',sub('.*gene_name "','',gtf[,9]))
d <- d[rownames(d) %in% trid,]
d <- rowsum(d,gn[match(rownames(d),trid)])

gtf <- fread('/home/whou10/scratch16/whou10/resource/gtf/grch38.gtf',data.table=F)
gtf <- gtf[gtf[,3]=='exon',]
gn <- sub('".*','',sub('.*gene_name "','',gtf[,9]))
gl <- tapply(gtf[,5]-gtf[,4],list(gn),sum)
glv <- as.vector(gl)
names(glv) <- names(gl)
d <- d/glv[rownames(d)]*1e3
d <- t(t(d)/colSums(d)*1e6)
d <- log2(d + 1)
saveRDS(d,file='/home/whou10/data/whou/metpred/data/hsc/proc/rna_matchedtrain.rds')
saveRDS(w,file='/home/whou10/data/whou/metpred/data/hsc/proc/wgbs_matchedtrain.rds')



