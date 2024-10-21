rm(list=ls())
ddir = '/home/whou10/data/whou/metpred/evaluate/encode/encode_cv/pred/ramp/'
rdir = '/home/whou10/data/whou/metpred/evaluate/encode/encode_cv/cgi/res/'

res = readRDS(paste0(rdir, 'mean_DNAm_of_cgi_that_overlap_with_promoter2e3Up.rds'))
res2 = readRDS(paste0(rdir, 'DNAm_sd_of_cgi_that_overlap_with_promoter2e3Up.rds'))
cgi.cor = readRDS(paste0(rdir, 'DNAm_cpg_cor_of_cgi_that_overlap_with_promoter2e3Up.rds'))

o = readRDS(paste0(rdir, 'o.df.rds'))
int = intersect(o[,2], rownames(res))
m=readRDS('/home/whou10/data/whou/metpred/evaluate/encode/cgi/res/mean_DNAm_of_cgi_that_overlap_with_promoter2e3Up.rds')
plot(rowMeans(m), pch = 20, cex = .3)
r = rowMeans(m)
o = o[match(rownames(m),o[,2]), ]
identical(o[,2], rownames(m))
ghigh =  o[r>0.9,1]
glow =  o[r<0.1,1]
tb = read.csv('/home/whou10/scratch16/whou10/resource/cancerGeneList.csv', sep = ',')
table(tb$Is.Oncogene, tb$Is.Tumor.Suppressor.Gene)
table(tb$Is.Oncogene)
gset = tb$Hugo.Symbol
names(gset) = gset
gset[tb$Is.Oncogene=='Yes' & tb$Is.Tumor.Suppressor.Gene=='Yes'] <- 'Onco.Suppressor'
gset[tb$Is.Oncogene=='Yes' & tb$Is.Tumor.Suppressor.Gene=='No'] <- 'Onco.nonSuppressor'
gset[tb$Is.Oncogene!='Yes'] <- 'nonOnco'
int = intersect(o[,1], names(gset))

rownames(m) = o[,1]
gset = gset[names(gset) %in% int]

gsetall <- rep('nonOnco', nrow(m))
names(gsetall) <- rownames(m)
gsetall[int] <- gset[int]
tab = table(gsetall)
for (i in 1:length(tab)){
  gsetall[gsetall==names(tab)[i]] <- paste0(names(tab)[i], '(', tab[i], ')')
}
  

pd = data.frame(gene = names(gsetall),
                Oncostatus = gsetall,
                dnam = rowMeans(m),
                stringsAsFactors = FALSE)

pdir <- '/home/whou10/data/whou/metpred/evaluate/encode/cgi/plot/'
library(ggplot2)
source('/home/whou10/scratch16/whou10/resource/myfunc/ggplot_theme.R')
pdf(paste0(pdir, 'mean_dnam_of_cgi_oncogene.pdf'), width = 5, height = 2)
ggplot(data = pd, aes(x = Oncostatus, y = dnam)) +
  geom_violin(aes(fill = Oncostatus))+
  theme_classic() +
  ylab('Mean DNAm of CGI') + 
  scale_fill_brewer(palette = 'RdYlBu') +
  theme(legend.position = 'none')
dev.off()

