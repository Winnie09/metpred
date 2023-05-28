## ENCODE
d.encode <- readRDS('/home/whou10/data/whou/metpred/data/encode_wgbs/processed/GRCh38_normal/rna.rds')
d.gtex <- readRDS('/home/whou10/data/whou/metpred/data/gtex_wgbs/processed/rna.rds') # 9
## /home/whou10/data/whou/metpred/data/v_20220624/combine 456 (Blueprint: 116, ENCODE: 96, GEO: 244) [exclude encode]
d.blueprint <- readRDS('/home/whou10/data/whou/metpred/data/v_20220624/combine/rna/hg38_blueprint.rds')
af <- list.files('/home/whou10/data/whou/metpred/data/v_20220624/combine/rna/', pattern = 'hg19_GSE')
d.tmp <- lapply(af, function(f){
  d.tmp <- readRDS(paste0('/home/whou10/data/whou/metpred/data/v_20220624/combine/rna/', f))
  colnames(d.tmp)
})
str(d.tmp)    
d.ceehrc <- readRDS('/home/whou10/data/whou/metpred/data/v_2021_notused/CEEHRC_wgbs/processed/rna.rds')

m.encode <- data.frame(DataSource = 'ENCODE', 
                       TissueName = colnames(d.encode),
                       stringsAsFactors = F)
m.gtex <- data.frame(DataSource = 'GTEx',
                     TissueName = colnames(d.gtex),
                     stringsAsFactors = F)
m.blueprint <- data.frame(DataSource = 'Blueprint',
                          TissueName = colnames(d.blueprint),
                          stringsAsFactors = F)
m.ceehrc <- data.frame(DataSource = 'CEEHRC',
                       TissueName = colnames(d.ceehrc),
                      stringsAsFactors = F)

m.geo <- data.frame(DataSource = 'Gene Expression Omnibus', 
                    TissueName = unlist(d.tmp),
                    stringsAsFactors = F)
str(m.geo)

meta <- rbind(m.encode,
              m.gtex,
              m.blueprint,
              m.ceehrc,
              m.geo)
saveRDS(meta, '/home/whou10/data/whou/metpred/data/processed_meta/processed_meta_colnames.rds')
write.csv(meta, '/home/whou10/data/whou/metpred/data/processed_meta/processed_meta_colnames.csv',
          row.names = F)

