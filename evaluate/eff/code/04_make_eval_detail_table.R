setwd('/home/whou10/data/whou/metpred/')
t = readRDS('./evaluate/eff/res/time.rds')
t = as.numeric(t[-1])
numberSample = seq(1e2, 1e3+200, 100)
names(t) = paste0('time_', numberSample,'(min)')
t_scaled = readRDS('./evaluate/eff/res/time_scaled.rds')
names(t_scaled) = paste0('time_',numberSample,'_scaled')
time = readRDS('./evaluate/eff/res/time.rds')

m = readRDS('./evaluate/eff/res/memory.rds')
m = m[-1]
m = as.numeric(sub('K','',m))/(1024^2)
names(m) = paste0('memory_', numberSample,'(GB)')

m_scaled = readRDS('./evaluate/eff/res/memory_scaled.rds')
names(m_scaled) = paste0('memory_', numberSample, '_scaled')
memory = readRDS('./evaluate/eff/res/memory.rds')
scalability = readRDS('./evaluate/eff/res/scalability.rds')

o = names(t)
df = data.frame(numberSample = numberSample, 
                t=t, t_scaled=t_scaled, time=time[-1], m=m, m_scaled=m_scaled, memory=memory[-1], scalability=scalability)

write.csv(df,'./evaluate/eff/res/efficiency_detail_scores.csv')

