proj <- readRDS('/home-4/zji4@jhu.edu/scratch/metpred/data/tcga/proc/combine/project.rds')
proj = gsub('TCGA-','',proj)

library(ggplot2)
pd = data.frame(proj = proj)
dir.create('/home-4/whou10@jhu.edu/scratch/Wenpin/metpred/model/lm/plot/stat/')
pdf('/home-4/whou10@jhu.edu/scratch/Wenpin/metpred/model/lm/plot/stat/study_abbreviation_barplot.pdf',width=8,height=3.5)
ggplot(pd, aes(x=proj)) + stat_count() + 
  geom_text(stat='count', aes(label=..count..), vjust=-1, size=3)+
  theme(axis.text.x = element_text(angle=90, hjust=1, size=10)) + xlab('Study abbreviation')
dev.off()

