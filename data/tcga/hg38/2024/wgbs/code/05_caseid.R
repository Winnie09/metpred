library(TCGAbiolinks)

pro <- TCGAbiolinks:::getGDCprojects()$project_id
pro <- grep('TCGA-',pro,value=T)

# query <- GDCquery(project = pro,
#                   legacy = TRUE,
#                   data.category = "DNA methylation",
#                   data.type = "Methylation percentage",
#                   experimental.strategy = "Bisulfite-Seq")
# res <- query$results[[1]]

query <- GDCquery(project = pro,
                  data.category = "DNA Methylation",
                  platform = c("Illumina Human Methylation 450"))
me <- query$results[[1]]

query <- GDCquery(project = pro,
                  data.category = "Transcriptome Profiling", 
                  data.type = "Gene Expression Quantification") 
 #                 workflow.type = "HTSeq - FPKM-UQ")
ge <- query$results[[1]]

me$patient <- me$cases.submitter_id
ge$patient <- ge$cases.submitter_id

me <- me[me$data_type=='Methylation Beta Value',]
saveRDS(me,'/home/whou10/data/whou/metpred/data/tcga/hg38/wgbs/meta/me.rds')
saveRDS(ge,'/home/whou10/data/whou/metpred/data/tcga/hg38/wgbs/meta/ge.rds')

