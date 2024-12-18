## gene list
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


## directory of predicted promoter DNAm profiles
ddir = '/home/whou10/data/whou/metpred/evaluate/mousepancreas/spatialpred/promoterDNAm/res/'
af = list.files(ddir)
af = setdiff(af,'df_of_cpg_gene.rds')
## check if the cell numbers match
ge = sapply(af, function(f){
  d = readRDS(paste0(ddir, f))
  # t(d[ag,])
  ag[ag%in% rownames(d)]
})


ge = do.call(rbind,sapply(af, function(f){
  d = readRDS(paste0(ddir, f))
  t(d[ag,])
}))
rownames(ge) <- sub('.*:','',rownames(ge))


## predicted promoter DNAm profiles have more cells:   
## A34_M_P_K4_D2.rds.2 (1 more) 
## A60_M_WT_LM_D2.rds.4 (1 more)
## A31_M_P_K4_D4.rds.5 (5 more)


## directory of results
rdir = '/home/whou10/data/whou/metpred/evaluate/mousepancreas/spatialpred/promoterDNAm/seuratobj/'
load('/home/whou10/data/whou/FeinbergLab/mousePancreas/visium/data/adm_log-normalized_with_sc_reference.rda')
int <- intersect(rownames(ge),rownames(adm.combined@meta.data))
ige <- matrix(NA,nrow=nrow(adm.combined@meta.data),ncol=ncol(ge),dimnames=list(row.names(adm.combined@meta.data),colnames(ge)))
ige[int,] <- ge[int,]
colnames(ige) <- paste0('DNAm_',colnames(ige))
adm.combined@meta.data <- data.frame(adm.combined@meta.data, ige)
saveRDS(adm.combined, paste0(rdir, 'seuratobj_20slides.rds'))


