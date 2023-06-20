## evaluate spatial pred
setwd('/home/whou10/data/whou/metpred/evaluate/mousepancreas/')
source('/home/whou10/scratch16/whou10/resource/startup.R')
dir.pf <- 'spatialpred/perf/pcc/'
dir.pf2 <- 'spatialpred/perf/mse/'
dir.p <- 'spatialpred/plot/perf/pcc/'
dir.p2 <- 'spatialpred/plot/perf/mse/'
test.s <- 'A46_F_P_K4_D4'
pred <- readRDS(paste0('spatialpred/res/pred_', test.s, '.rds'))
genecpg <- readRDS('spatialpred/genecpg/gene_cpg_leaveout1k.rds')
predagg <- aggregate(pred, list(genecpg[match(rownames(pred), genecpg[,1]), 2]), 'mean')
rownames(predagg) <- predagg[,1]
predagg <- as.matrix(predagg[, -1])
expr <- readRDS('/home/whou10/data/whou/metpred/data/mousepancreas/spatial/normexpr.rds')
gs <- expr[rownames(predagg), colnames(predagg)]
gs <- as.matrix(gs)
summary(rowsds(gs))

## evaluation
cc <- corfunc(predagg, gs)
saveRDS(cc, paste0(dir.pf, test.s, '.rds'))

cc.acrossgene <- corfunc(t(predagg), t(gs))
saveRDS(cc.acrossgene, paste0(dir.pf, test.s, 'acrossgene_within_spot.rds'))

mse <- sapply(1:nrow(predagg), function(i){
  tmpv <- predagg[i, ] - gs[i, ]
  mean(tmpv * tmpv)
})
names(mse) <- rownames(predagg)
saveRDS(mse, paste0(dir.pf2, test.s, '.rds'))

### plot
rs <- rowsds(gs)
rs2 <- cut(rs, seq(0, 1, 0.05))
names(rs2) <- names(rs)
cc <- cc[names(rs2)]

pdf(paste0(dir.p, test.s ,'_acrossspot_pcc.pdf'), width=2.5,height=2.5)
ggplot(data.frame(sd=rs2,cor=cc),aes(y=cor,x=sd,fill=rs2)) + 
  geom_violin(alpha=0.2, scale = 'area') + 
  geom_boxplot(alpha=0.3, width = 0.2) + 
  theme_classic() + 
  scale_fill_manual(values=rainbow(length(unique(cc)))) + 
  xlab('Testing measured value sd') + 
  ylab('across-spot PCC') + 
  theme(legend.position = 'none', axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

pdf(paste0(dir.p2, test.s,'_acrossspot_mse.pdf'),width=2.5,height=2.5)
ggplot(data.frame(sd=rs2,mse=mse),aes(y=mse,x=sd,fill=rs2)) + 
  geom_violin(alpha=0.2, scale = 'area') + 
  geom_boxplot(alpha=0.3, width = 0.2) + 
  theme_classic() + 
  scale_fill_manual(values=rainbow(length(unique(cc)))) + 
  xlab('Testing measured value sd') + 
  ylab('across-spot MSE') + 
  theme(legend.position = 'none', axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

### plot
rs.acrossgene <- rowsds(t(gs))
rs.acrossgene2 <- cut(rs.acrossgene, seq(0, 1, 0.05))
names(rs.acrossgene2) <- names(rs.acrossgene)
cc.acrossgene2 <- cc.acrossgene[names(rs.acrossgene2)]

pdf(paste0(dir.p, test.s ,'_acrossgene_pcc.pdf'),width=2.5,height=2.5)
ggplot(data.frame(sd=rs.acrossgene2,cor=cc.acrossgene),aes(y=cor,x=sd,fill=rs.acrossgene2)) + 
  geom_violin(alpha=0.2, scale = 'area') + 
  geom_boxplot(alpha=0.3, width = 0.2) + 
  theme_classic() + 
  scale_fill_manual(values=rainbow(length(unique(cc.acrossgene2)))) + 
  xlab('Testing measured value sd') + 
  ylab('across-gene PCC') + 
  theme(legend.position = 'none', axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()


