#############
rdir = '/home/whou10/data/whou/metpred/evaluate/mousepancreas/spatialpred/promoterDNAm/seuratobj/'
pdir = '/home/whou10/data/whou/metpred/evaluate/mousepancreas/spatialpred/promoterDNAm/plot/'
## gene list
adm.combined = readRDS(paste0(rdir, 'seuratobj_20slides.rds'))
## plot seurat object
library(Seurat)

ag = c('Ins1',
       'Ins2',
       'Gcg',
       'Ppy',
       'Cpa1',
       'Mafa',
       'Klf4',
       'Klf5',
       'Hnf1b',
       'Hes1',
       'Prss1',
       'Pdx1',
       'Sst',
       'Epcam',
       'Krt19',
       'Sox9'
)

for (g in ag) {
  pdf(paste0(pdir,g,'.pdf'),width=15,height=15)
  plot(SpatialFeaturePlot(adm.combined,feature=paste0('DNAm_', g),ncol=5))
  dev.off()
}


