## global 
rdir <- '/home/whou10/data/whou/metpred/evaluate/encode/encode_cv/cgi/res/'
pdir <- '/home/whou10/data/whou/metpred/evaluate/encode/encode_cv/cgi/plot/'

## read in CpG island (CGI) location information
tb = read.table('/home/whou10/data/whou/metpred/data/cpg_island/hgTables.txt', stringsAsFactors = FALSE)
library(GenomicRanges)
gr <- GRanges(seqnames=tb[,2],IRanges(start=tb[,3],end=tb[,4]))

## read in predicted CpG location information of cancer samples
pred <- readRDS('/home/whou10/data/whou/metpred/evaluate/encode/encode_cv/pred/ramp/1.rds')
                # /home-4/whou10@jhu.edu/scratch/Wenpin/metpred/predict/res/predicted_DNAm_on_testset_cancer.rds')
rn <- rownames(pred)
start <- as.numeric(sub('.*:', '', rn))
cpg <- GRanges(seqnames= sub(':.*', '', rn), IRanges(start=start,end=start))

## find the overlap of CGI and CpG 
o <- as.data.frame(findOverlaps(gr,cpg))
table(table(o[,1]))


pdf(paste0(pdir, 'numberCpG_in_numberCpG_hist.pdf'), width = 6, height = 4)
hist(table(o[,1]), breaks = 50, col = 'grey', xlab = 'number of CpG', ylab = 'number of CGI (hg38)', main = paste0('Total: ', length(gr), ' CGI'))
dev.off()

## calculate the PCC  among the predicted profiles of the CpG within a same CGI
cgi.mean <- sapply(unique(o[,1]), function(i){
  tmp = o[o[,1] == i, , drop = FALSE]
  c = as.vector(cor(t(pred[tmp[,2], ])))
  mean(c[c!=1])
}) 
cgi.mean = cgi.mean[!is.na(cgi.mean)]

## permute the CGI-CpG relationship and recalculate the PCC among the predicted profiles of the CpG within a CGI
o.pm <- o
set.seed(12345)
o.pm[,2] <- sample(o[,2])
cgipm.mean <- sapply(unique(o[,1]), function(i){
  tmp = o.pm[o.pm[,1] == i, , drop = FALSE]
  c = as.vector(cor(t(pred[tmp[,2], ])))
  mean(c[c!=1])
}) 
cgipm.mean = cgipm.mean[!is.na(cgipm.mean)]


## read in encode samples measured DNA methylation (true values)

me <-
  readRDS(
    '/home/whou10/data/whou/metpred/data/encode_wgbs/processed/GRCh38_normal/wgbs.rds'
  )
str(me)
str(pred)
true = me[, colnames(pred)]

str(true)
# tb = read.csv('/home-4/whou10@jhu.edu/scratch/Wenpin/metpred/doc/ct_tissue_info_for_matched_wgbs_and_RNAseq.csv', header = TRUE, as.is = TRUE)
# m.bak <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/metpred/final/nonawgbs_hg38.rds')
# id = tb[tb[,3] == TRUE, 1] ## cancer
# str(id)
# m = m.bak[, id]
# testid <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/metpred/predict/res/testid_cancer.rds')
# true <- m[,testid]

## calculate the true DNAm  PCC among the CpG within a same CGI
cgi.mean.true <- sapply(unique(o[,1]), function(i){
  tmp = o[o[,1] == i, , drop = FALSE]
  c = as.vector(cor(t(true[tmp[,2], ])))
  mean(c[c!=1])
}) 
cgi.mean.true = cgi.mean.true[!is.na(cgi.mean.true)]




## compare and plot the CpG PCC within a same CGI using predicted, measured values in matched CpG-CGI relationship, and permuted CpG-CGI relationship
pd = data.frame(type = c(rep('pred.mean', length(cgi.mean)), rep('pmpred.mean', length(cgipm.mean)), rep('true.mean', length(cgi.mean.true))),
                value = c(cgi.mean, cgipm.mean, cgi.mean.true))
saveRDS(pd, paste0(rdir, 'pd.rds'))

pd = readRDS(paste0(rdir, 'pd.rds'))
pd[,1] = factor(pd[,1])
levels(pd[,1]) = c('Permute', 'Ramp', 'Measured')
colv = readRDS('/home/whou10/data/whou/metpred/setting/method_color.rds')
library(ggplot2)
library(RColorBrewer)
pdf(paste0(pdir, 'cgi_mean_violin.pdf'), width = 4, height = 2.7)
ggplot(data = pd, aes(x = type, y = value, fill = type)) + 
  geom_jitter(alpha = 0.1, stroke = 0, size = 0.1, width = 0.1) + 
  geom_violin(scale = 'width',width = 0.8, alpha = 0.3) + 
  geom_boxplot(outlier.shape = NA, width = 0.1) + 
  scale_fill_manual(values = colv) + 
  theme_classic() +
  theme(legend.position = 'none', 
        axis.text = element_text(color = 'black', size = 10), 
        axis.title = element_text(color = 'black', size = 10))+
  ylab('Mean DNAme within CpG islands') + 
  xlab('Models')
dev.off()

## find a CGI with more than 30 CpG as an example
v <- sapply(unique(o[,1]), function(i){
  tmp = o[o[,1] == i, , drop = FALSE]
  if (nrow(tmp) > 30){
    return(i)
  } else {
    return(NA)
  }
}) 
v2 <- v[!is.na(v)]

## plot the heatmap of the values (rows: CpG, columns: samples)
i = v2[1]
tmp = o[o[,1] == i, , drop = FALSE]
tmp.pm = o.pm[o.pm[,1] == i, , drop = FALSE]

hm.pd = rbind(pred[tmp[,2],],
              pred[tmp.pm[,2],],
              true[tmp[,2],])
rownames(hm.pd) <- paste0(rownames(hm.pd), ';', seq(1, nrow(hm.pd)))
rowann = data.frame(type = rep(c('pred', 'pred.pm', 'true'), each = nrow(tmp)))
rownames(rowann) = rownames(hm.pd)
rowann.col = c('orange', 'blue', 'purple')
names(rowann.col) = unique(rowann)

library(pheatmap)

library(grid)
pdf(paste0(pdir, 'example_CGI_predvalues_pmpredvalues_truevalues_hm.pdf'), width=5, height = 5)
## Create the heatmap:
setHook("grid.newpage", function() pushViewport(viewport(x=1,y=1,width=0.9, height=0.9, name="vp", just=c("right","top"))), action="prepend")
pheatmap(hm.pd, scale = 'none', show_rownames = F, show_colnames = F, cluster_rows = FALSE, annotation_row = rowann, border_color = NA, main = paste0(seqnames(gr)[i], ':', start(gr)[i], '-', end(gr)[i])) 
setHook("grid.newpage", NULL, "replace")
grid.text("Samples", y=0.01, gp=gpar(fontsize=10))
grid.text("CpGs", x=0, rot=90, gp=gpar(fontsize=10))
dev.off()


