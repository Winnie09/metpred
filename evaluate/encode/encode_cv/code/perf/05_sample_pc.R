### ENCODE predicted and measured sample PCA

## load data
#p <- readRDS('/home/whou10/data/whou/metpred/evaluate/encode_gtex/pred/ramp.rds')
#d <- readRDS('/home/whou10/data/whou/metpred/data/encode_wgbs/processed/GRCh38_normal/wgbs.rds')

cvfold = commandArgs(trailingOnly = T)[[1]][1]
print(cvfold)
p <- readRDS(paste0('/home/whou10/data/whou/metpred/evaluate/encode_cv/pred/ramp/', cvfold, '.rds'))
# > str(p)
# num [1:28301739, 1:10] 0.817 0.815 0.804 0.799 0.796 ...
# - attr(*, "dimnames")=List of 2
# ..$ : chr [1:28301739] "chr1:10469" "chr1:10471" "chr1:10484" "chr1:10489" ...
# ..$ : chr [1:10] "Homo sapiens aorta tissue male adult (34 years)" "Homo sapiens urinary bladder tissue male child (3 years)" "Homo sapiens B cell male adult (37 years)" "Homo sapiens CD14-positive monocyte male adult (37 years)" ...

d <- readRDS('/home/whou10/data/whou/metpred/data/encode_wgbs/processed/GRCh38_normal/wgbs.rds')
# > str(d)
# num [1:28301739, 1:95] 0.827 0.824 0.805 0.798 0.793 ...
# - attr(*, "dimnames")=List of 2
# ..$ : chr [1:28301739] "chr1:10469" "chr1:10471" "chr1:10484" "chr1:10489" ...
# ..$ : chr [1:95] "Homo sapiens aorta tissue male adult (34 years)" "Homo sapiens urinary bladder tissue male child (3 years)" "Homo sapiens B cell male adult (37 years)" "Homo sapiens CD14-positive monocyte male adult (37 years)" ...

str(p)
str(d)
colnames(p) = paste0('Pred:', colnames(p))

d <- cbind(d,p)

## simplify the cell type /tissue names
ct <- rep('other',ncol(d))
ct[colnames(d) %in% c(sub('Pred:','',colnames(p)), colnames(p))] <- colnames(d)[colnames(d) %in% c(sub('Pred:','',colnames(p)), colnames(p))]
ct = sub('Homo sapiens ', '', sub('Pred:','',sub(' tissue.*','', ct)))

# ct[grep('Esophagus.?muscularis',colnames(d),ignore.case = T)] <- 'esophagus\nmuscularis'
ct[grep('heart',colnames(d),ignore.case = T)] <- 'heart'
ct[grep('lung',colnames(d),ignore.case = T)] <- 'lung'
# ct[grep('skeletal',colnames(d),ignore.case = T)] <- 'skeletal\nmuscle'
ct[grep('skin',colnames(d),ignore.case = T)] <- 'skin'
ct[grep('thyroid',colnames(d),ignore.case = T)] <- 'thyroid'
# ct[grep('urinary bladder',colnames(d),ignore.case = T)] <- 'urinary\nbladder'
# ct[grep('CD14-positive monocyte',colnames(d),ignore.case = T)] <- 'CD14+ monocyte'
# ct[grep('common myeloid progenitor, CD34-positive',colnames(d),ignore.case = T)] <- 'CD34+ myeloid\nprogenitor'
ct[grep('esophagus squamous epithelium',colnames(d),ignore.case = T)] <- 'ESE'
ct[grep('gastroesophageal sphincter',colnames(d),ignore.case = T)] <- 'gastroesophageal\nsphincter'
ct[grep('esophagus tissue',colnames(d),ignore.case = T)] <- 'esophagus'
ct[grep('spleen',colnames(d),ignore.case = T)] <- 'spleen'
ct[grep('adipose',colnames(d),ignore.case = T)] <- 'adipose'
ct[grep('adrenal gland',colnames(d),ignore.case = T)] <- 'adrenal\ngland'
ct[grep('aorta',colnames(d),ignore.case = T)] <- 'aorta'
ct[grep('T-cell',colnames(d),ignore.case = T)] <- 'T cell'
ct[grep('B cell',colnames(d),ignore.case = T)] <- 'B cell'
ct[grep('GM',colnames(d),ignore.case = T)] <- 'other'

sort(ct)
length(unique(ct))


## only select large tissue types
d.bak = d
d <- d[,ct!='other']
ct <- ct[ct!='other']


## PCA
rowsds <- function(data,cm) {
  sqrt((rowMeans(data*data) - cm^2) / (ncol(data) - 1) * ncol(data))
}
v = rowsds(d,rowMeans(d))
summary(v)
d3 <- d[v > 0.1,]
d = d3
pr <- prcomp(t(d),scale. = T)$x
saveRDS(pr, paste0('/home/whou10/data/whou/metpred/evaluate/encode_cv/plot/cv', cvfold,'_pr.rds'))

db <- ifelse(grepl('^Pred',rownames(pr)),'Predicted','Measured')

## plot
library(ggplot2)
library(ggrepel)
library(RColorBrewer)
label = ct
label[duplicated(ct)] <- ''
pd <- data.frame(PC1=pr[,1],PC2=pr[,2],ct=ct,db=db,label=label,stringsAsFactors = F)
saveRDS(pd, paste0('/home/whou10/data/whou/metpred/evaluate/encode_cv/plot/cv', cvfold,'_pd.rds'))

colv = readRDS('/home/whou10/scratch16/whou10/resource/color28.rds')
colv = colv[1:length(unique(ct))]
names(colv) = sort(unique(ct))
pdf(paste0('/home/whou10/data/whou/metpred/evaluate/encode_cv/plot/cv', cvfold,'_pc.pdf'),width=4.2,height=4.2)
ggplot(pd,aes(x=PC1,y=PC2,col=ct,shape=db,label=label)) + 
  geom_point(size=2) + 
  geom_text_repel(max.overlaps = 10) + 
  theme_classic() + 
  theme(legend.position = 'none',legend.title = element_blank())+
  xlab(paste0('Principal component 1 (', round(var(pr[,1])/sum(apply(pr,2,var)) * 100,2), '%)'))+
  ylab(paste0('Principal component 2 (', round(var(pr[,2])/sum(apply(pr,2,var)) * 100,2), '%)'))+
  scale_color_manual(values = colv)
dev.off()


