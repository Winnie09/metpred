## 20230530 mouse pancreas, 
## visium's deconvoluted celltypes (duct, acinar, adm) distribution for each visium sample
rdir <- '/home/whou10/data/whou/metpred/evaluate/mousepancreas/goldstandard/res/'
pdir <- '/home/whou10/data/whou/metpred/evaluate/mousepancreas/goldstandard/plot/'

w <- readRDS('/home/whou10/data/whou/metpred/data/mousepancreas/spatial/weight.rds')
str(w)
alls <- sub(':.*', '', rownames(w))
length(unique(alls))

## each Visium sample, the WGBS cell types
allf <- list.files('/home/whou10/data/whou/FeinbergLab/mousePancreas/wgbs/data/bamfile')
v <- sub('_aci.*','', allf)
v <- sub('_duc.*','', v)
v <- sub('_ADM.*','', v)
v <- sub('_unr.*','', v)
ctlist <- list()
for (i in v){
  tmp <- allf[v == i]
  ctlist[[i]] <- sapply(tmp, function(j){
    sub('\\..*', '', sub(paste0(i, '_'), '', j))
  }, USE.NAMES = F)  
}
saveRDS(ctlist, paste0(rdir, 'wgbs_celltypes_of_each_visium_sample.rds'))

## the deconvolute expr weight of WGBS cell types
wlist <- list()
for (s in unique(alls)){
 print(s)
 ct.s <- tolower(ctlist[[s]])
 if ("adm" %in% ct.s) ct.s <- c(ct.s, 'severe_adm')
 if ('unrecovered_adm' %in% ct.s)  ct.s <- setdiff(ct.s, 'unrecovered_adm')
 tmp <- w[alls == s, ct.s, drop = F] 
 
 if ('adm' %in% ct.s) {
   tmp[,'adm'] <- tmp[,'adm'] + tmp[,'severe_adm']
   tmp <- tmp[, !colnames(tmp) %in% c('severe_adm')]
 }
 wlist[[s]] <- tmp
}



## plot the deconvolute expr weight of WGBS cell types
pd.all <- lapply(1:length(wlist), function(i){
  m = as.matrix(wlist[[i]])
  str(m)
  pd <- reshape2::melt(m)
  colnames(pd) <- c('spot', 'celltype', 'weight')
  pd$sample <- sub(':.*', '', pd[,1])
  pd
})
str(pd.all)  
pd.all <- do.call(rbind, pd.all)


library(ggplot2)
library(RColorBrewer)


# Draw two histograms in same plot
pdf(paste0(pdir, 'deconvoluted_weight_hist.pdf'), height = 4, width = 6)
ggplot(pd.all, aes(x = weight, fill = celltype, color = celltype)) +            
  geom_histogram(alpha = 0.3, size = 0.2, position = "identity", bins = 20) + 
  geom_freqpoly(alpha = 1, bins = 20) + 
  scale_color_brewer(palette = 'Set2') +
  scale_fill_brewer(palette = 'Pastel2') + 
  xlab('deconvoluted weight') + 
  ylab('number of spots')  + 
  facet_wrap(~sample)
dev.off()


## only calculate the sum of the ct weights that have wgbs data
## select spots with sum > 0.9
### first, calculate the number of spots
numberOfSpot <- sapply(seq(0.8, 0.95, 0.05), function(cutoff){
  num.spot <- sapply(names(wlist), function(i){
    tmp <- wlist[[i]]
    sum(rowSums(tmp) > cutoff)
  })
})
str(numberOfSpot)  
colnames(numberOfSpot) <- paste0('weightSum>', seq(0.8, 0.95, 0.05))
write.csv(numberOfSpot, paste0(rdir, 'numberOfSpot_by_cutoffs.csv'))

### select spots
tmp <- sapply(seq(0.8, 0.95, 0.05), function(cutoff){
  spot.list <- lapply(names(wlist), function(i){
    tmp <- wlist[[i]]
    v <- rownames(tmp)[rowSums(tmp) > 0.95]
  })
  names(spot.list) <- names(wlist)
  saveRDS(spot.list, paste0(rdir, 'weightSum_cutoff_', cutoff, '_spotList.rds'))
  return(0)
})

