## ==================================================
## use TCGAbiolinks: retrieve DNAm samples. They are all 450 array data
## ==================================================
library(data.table)
library('TCGAbiolinks')
query <- GDCquery(
    project = c("TCGA-GBM", "TCGA-LGG"),
    data.category = "DNA Methylation",
    # platform = c("Illumina Human Methylation 450"),
    # sample.type = "Recurrent Tumor"
)
dt = data.table(
    getResults(query), 
    filter = 'top',
    options = list(scrollX = TRUE, keys = TRUE, pageLength = 5), 
    rownames = FALSE
)

table(dt$platform)

 # Illumina Human Methylation 27 Illumina Human Methylation 450 
 #                           885                           2067 
summary(dt$file_size)
  #  Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
  # 719548   766748  8095273  7059314  8095341 13214929 
 ## ==================================================
