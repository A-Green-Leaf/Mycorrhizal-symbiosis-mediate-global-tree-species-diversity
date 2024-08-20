library(MASS)
library(spdep)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(ggpp)
library(car)

global.pro <- readRDS('data/global.pro.rds')

source('script/functions.R')

summary(global.pro)
table(global.pro$Biome)

# Biomes: GLM model -------------------------------------------------------
col.use <- c(
  'aquamarine',
  'gray','gray1','gray20',
  'navajowhite','orange','orange3',
  'orangered','orangered2','orangered4',
  'aquamarine3')

data.biome <- global.pro %>% 
  dplyr::select(SR, EcM.p.stem, EcM.p.stem2, bio1,AI,elev,Slope,sg.clay.0_200,sg.cec.0_200,DBH,Area_log, Lat_abs,Biome) %>% na.omit()

biome <- table(data.biome$Biome)[which(table(data.biome$Biome)>=100)] %>% names()
biome.list <- data.frame(Biome=biome, col.use=col.use)

biome.name <- gsub(' ','_',biome)
biome.name <- gsub('/','_',biome.name)
biome.name <- gsub('&','',biome.name)
  
  
# model
biome.mod <- list()
for (i in 1:length(biome)) {
  data.i <- subset(data.biome, Biome==biome[i])
  biome.mod[[i]] <- glm.biome(data.i, biome[i])
  print(paste0(i, '_ overdispersion is ', biome.mod[[i]][[4]]))
}
names(biome.mod) <- biome.name



# Biomes: Figure --------------------------------------------------------------


# predict data
pred.biome <- data.frame()
for (i in 1:length(biome)) {
  out.i <- biome.mod[[i]]$newdata
  pred.biome <- rbind(out.i, pred.biome)
}

# coef data
coef.biome <- data.frame()
for (i in 1:length(biome)) {
  out.i <- biome.mod[[i]]$coef
  coef.biome <- rbind(out.i, coef.biome)
}


# plot lines

# biome legend
biome.legend = ggplot()+
  geom_point(aes(x=rep(1,11), y=seq(2,10, length=11)), shape=21,fill=col.use)+
  annotate('text',x=rep(1.5,11), y=seq(2,10, length=11), label=biome, hjust = 0,size=2)+
  xlim(0,11)+ylim(0,11)+
  theme_void()

# 1. biome predicted lines
a.biome=ggplot(pred.biome, aes(EcM.p.stem, pred, col=name))+
  geom_line(linewidth=0.6)+scale_y_log10()+
  scale_color_manual(values = col.use)+
  theme_classic()+xlab('Proportion of EcM tree abundance')+ylab('Species richness')+
  theme(legend.position = 'none',#plot.margin = ggplot2::margin(0.5, 0, 0, -0.2, "cm"),axis.title.y = element_text(vjust=-2.3),
     plot.title = element_text(face = 'bold'))+ggtitle('a')

# 2. Linear slopes vs. latitude and AI
with(subset(coef.biome, variable=='EcM.p.stem'), cor.test(Lat_abs, Estimate))
with(subset(coef.biome, variable=='EcM.p.stem'), cor.test(AI, Estimate))

b.biome=ggplot(subset(coef.biome, variable=='EcM.p.stem'), aes(Lat_abs, Estimate))+
  geom_smooth(method = glm, col='black',se=F)+
  geom_point(shape=21, size=3,aes(fill=name))+
  scale_fill_manual(values = col.use)+
  theme_classic()+xlab('Absolute latitude')+ylab('Linear slope')+
  theme(legend.position = 'none', plot.title = element_text(face = 'bold'))+ggtitle('b')+
  geom_text_npc(aes(npcx =0.95, npcy =0.15, label='italic(r)==0.797'),parse=TRUE)+
  geom_text_npc(aes(npcx =0.95, npcy =0.05, label='italic(p)==0.003'),parse=TRUE)

b1.biome=ggplot(subset(coef.biome, variable=='EcM.p.stem'), aes(AI, Estimate))+
  geom_smooth(method = glm, col='black',se=F)+
  geom_point(shape=21, size=3,aes(fill=name))+
  scale_fill_manual(values = col.use)+
  theme_classic()+xlab('Absolute latitude')+ylab('Linear slope')+
  theme(legend.position = 'none', plot.title = element_text(face = 'bold'))+ggtitle('b')+
  geom_text_npc(aes(npcx =0.95, npcy =0.15, label='italic(r)==0.797'),parse=TRUE)+
  geom_text_npc(aes(npcx =0.95, npcy =0.05, label='italic(p)==0.003'),parse=TRUE)

# 3. Quadratic slopes vs. latitude and AI
with(subset(coef.biome, variable=='EcM.p.stem2'), cor.test(Lat_abs, Estimate))
with(subset(coef.biome, variable=='EcM.p.stem2'), cor.test(AI, Estimate))

c.biome=ggplot(subset(coef.biome, variable=='EcM.p.stem2'), aes(Lat_abs, Estimate))+
  geom_smooth(method = glm, col='black',se=F)+
  geom_point(shape=21, size=3,aes(fill=name))+
  scale_fill_manual(values = col.use)+
  theme_classic()+xlab('Absolute latitude')+ylab('Quadratic slope')+
  theme(legend.position = 'none', plot.title = element_text(face = 'bold'))+ggtitle('c')+
  geom_text_npc(aes(npcx =0.95, npcy =0.15, label='italic(r)==-0.707'),parse=TRUE)+
  geom_text_npc(aes(npcx =0.95, npcy =0.05, label='italic(p)==0.015'),parse=TRUE)

c1.biome=ggplot(subset(coef.biome, variable=='EcM.p.stem2'), aes(AI, Estimate))+
  geom_smooth(method = glm, col='black',se=F)+
  geom_point(shape=21, size=3,aes(fill=name))+
  scale_fill_manual(values = col.use)+
  theme_classic()+xlab('Absolute latitude')+ylab('Quadratic slope')+
  theme(legend.position = 'none', plot.title = element_text(face = 'bold'))+ggtitle('c')+
  geom_text_npc(aes(npcx =0.95, npcy =0.15, label='italic(r)==-0.707'),parse=TRUE)+
  geom_text_npc(aes(npcx =0.95, npcy =0.05, label='italic(p)==0.015'),parse=TRUE)

# 4. Linear + Quadratic vs. latitude and AI
x.biome=subset(coef.biome, variable=='EcM.p.stem') %>%
  mutate(
    linear = Estimate,
    quadratic = subset(coef.biome, variable=='EcM.p.stem2')$Estimate,
    Diff =  linear+quadratic) %>%
  mutate(
    hypo='positive unimodal',
    hypo=ifelse(Diff>0&quadratic>=0, 'positive',hypo),
    hypo=ifelse(Diff<0, 'Negtive unimodal',hypo),
    hypo=ifelse(Diff<0&linear<=0, 'Negtive',hypo)
  )

with(x.biome, cor.test(Lat_abs, Diff))
with(x.biome, cor.test(AI, Diff))

factor(x.biome$hypo)
diff.biome <- ggplot(x.biome, aes(Lat_abs, Diff))+
  geom_smooth(method = glm, col='black',se=F)+
  geom_point(shape=21, size=3,aes(fill=hypo))+
  scale_fill_manual(values = c('#F8766D', '#C77CFF', '#7CAE00'))+
  theme_classic()+xlab('Absolute latitude')+ylab('Linear+Quadratic slopes')+
  theme(legend.position = 'none', legend.title = element_blank(), legend.background = element_blank(), plot.title = element_text(face = 'bold'))+ggtitle('b')+
  geom_hline(yintercept = 0, lty=2, col='red')+
  geom_text_npc(aes(npcx =0.95, npcy =0.15, label='italic(r)==0.874'),parse=TRUE)+
  geom_text_npc(aes(npcx =0.95, npcy =0.05, label='italic(p)<0.001'),parse=TRUE)+
  geom_point(data=data.frame(x=rep(12,4), y=seq(2.5,4, length=4)),aes(x=x, y=y), shape=21,fill=c('#F8766D', '#C77CFF', '#7CAE00', '#00BFC4'),size=2)+
  annotate('text',x=rep(14,4), y=seq(2.5,4, length=4), label=c('Negative','Negative unimodal','Posive unimodal','Positive'), hjust = 0,size=2)

diff1.biome <- ggplot(x.biome, aes(AI, Diff))+
  geom_smooth(method = glm, col='black',se=F)+
  geom_point(shape=21, size=3,aes(fill=hypo))+
  scale_fill_manual(values = c('#F8766D', '#C77CFF', '#7CAE00'))+
  theme_classic()+xlab('Aridity index')+ylab('Linear+Quadratic slopes')+
  theme(legend.position = 'none', legend.title = element_blank(), legend.background = element_blank(), plot.title = element_text(face = 'bold'))+ggtitle('c')+
  geom_hline(yintercept = 0, lty=2, col='red')+
  geom_text_npc(aes(npcx =0.95, npcy =0.15, label='italic(r)==-0.099'),parse=TRUE)+
  geom_text_npc(aes(npcx =0.95, npcy =0.05, label='italic(p)==0.772'),parse=TRUE)





# GLM across ecoregions ---------------------------------------------------

data.ecoregion <- global.pro %>% 
  dplyr::select(SR, EcM.p.stem, EcM.p.stem2, bio1,AI,elev,Slope,sg.clay.0_200,sg.cec.0_200,DBH,Area_log, Lat_abs,Ecoregion,Biome) %>% 
  na.omit() %>%
  mutate(
    bio1=bio1+rnorm(length(bio1), 0, 0.00001),
    AI=AI+rnorm(length(AI), 0, 0.00001),
    elev=elev+rnorm(length(elev), 0, 0.00001),
    Slope=Slope+rnorm(length(Slope), 0, 0.00001),
    sg.clay.0_200=sg.clay.0_200+rnorm(length(sg.clay.0_200), 0, 0.00001),
    sg.cec.0_200=sg.cec.0_200+rnorm(length(sg.cec.0_200), 0, 0.00001),
    DBH=DBH+rnorm(length(DBH), 0, 0.00001),
    Area_log=Area_log+rnorm(length(Area_log), 0, 0.00001)
  )

# select ecoregion
ecoregion.tab <- as.data.frame(table(data.ecoregion$Ecoregion))
colnames(ecoregion.tab) <- c('Ecoregion','N')
ecoregion.list <- data.ecoregion %>% 
  left_join(ecoregion.tab, by=c('Ecoregion')) %>%
  group_by(Biome, Ecoregion) %>%
  summarise(N=mean(N, na.rm=T), EcMrange=max(EcM.p.stem)-min(EcM.p.stem)) %>%
  filter(N>=100 & EcMrange>=0.5) %>%
  left_join(biome.list, by='Biome')

ecoregion <- ecoregion.list$Ecoregion

ecoregion.name <- gsub(' ','_',ecoregion)
ecoregion.name <- gsub('/','_',ecoregion.name)
ecoregion.name <- gsub('&','',ecoregion.name)


# model
ecoregion.mod <- list()
for (i in 1:length(ecoregion)) {
  data.i <- subset(data.ecoregion, Ecoregion==ecoregion[i])
  ecoregion.mod[[i]] <- glm.ecoregion(data.i, ecoregion[i])
  print(paste0(i, '_ overdispersion is ', ecoregion.mod[[i]][[4]]))
}
names(ecoregion.mod) <- ecoregion.name



# Plot ecoregions --------------------------------------------------------------


# predict data
pred.ecoregion <- data.frame()
for (i in 1:length(ecoregion)) {
  out.i <- ecoregion.mod[[i]]$newdata
  pred.ecoregion <- rbind(out.i, pred.ecoregion)
}
pred.ecoregion <- merge(pred.ecoregion, unique(global.pro[,c('Ecoregion','Biome')]), by.x='name', by.y='Ecoregion', all.x=T, sort=F)
names(table(pred.ecoregion$Biome))

# coef data
coef.ecoregion <- data.frame()
for (i in 1:length(ecoregion)) {
  out.i <- ecoregion.mod[[i]]$coef
  coef.ecoregion <- rbind(out.i, coef.ecoregion)
}
coef.ecoregion <- merge(coef.ecoregion, unique(global.pro[,c('Ecoregion','Biome')]), by.x='name', by.y='Ecoregion', all.x=T, sort=F)
names(table(coef.ecoregion$Biome))

col.ecoregion <- c(
  'aquamarine',
  'gray','gray1','gray10','gray20',
  'navajowhite','orange','orange3',
  'red','orangered','orangered2','orangered4',
  'aquamarine3'
)

col.ecoregion <- c(
  'aquamarine',
  'gray','gray1','gray20',
  'navajowhite','orange','orange3',
  'orangered','orangered2','orangered4',
  'aquamarine3'
)
# plot lines
# plot ecoregion legend
ecoregion.legend=ggplot()+
  geom_point(aes(x=rep(1,11), y=seq(2,10, length=11)), shape=21,fill=col.use)+
  annotate('text',x=rep(1.5,11), y=seq(2,10, length=11), label=biome, hjust = 0,size=2)+
  xlim(0,11)+ylim(0,11)+
  theme_void()

# predict lines
a.ecoregion=ggplot(pred.ecoregion, aes(EcM.p.stem, pred, group=name, col=Biome))+
  geom_line(linewidth=0.3,alpha=0.5)+scale_y_log10()+
  scale_color_manual(values = col.ecoregion)+
  theme_classic()+xlab('Proportion of EcM tree abundance')+ylab('Species richness')+
  theme(legend.position = 'none', plot.title = element_text(face = 'bold'))+ggtitle('a')

# 1. Linear slopes vs. Latitude and AI
with(subset(coef.ecoregion, variable=='EcM.p.stem'), cor.test(Lat_abs, Estimate))
with(subset(coef.ecoregion, variable=='EcM.p.stem'), cor.test(AI, Estimate))

b.ecoregion=ggplot(subset(coef.ecoregion, variable=='EcM.p.stem'), aes(Lat_abs, Estimate))+
  geom_smooth(method = glm, col='black',se=F)+
  geom_point(shape=21, size=3,aes(fill=Biome))+
  scale_fill_manual(values = col.ecoregion)+
  theme_classic()+xlab('Absolute latitude')+ylab('Linear slope')+
  theme(legend.position = 'none')+
  geom_text_npc(aes(npcx =0.05, npcy =0.75, label='italic(r)==0.482'),parse=TRUE)+
  geom_text_npc(aes(npcx =0.05, npcy =0.65, label='italic(p)<0.001'),parse=TRUE)

b1.ecoregion=ggplot(subset(coef.ecoregion, variable=='EcM.p.stem'), aes(AI, Estimate))+
  geom_smooth(method = glm, col='black',se=F)+
  geom_point(shape=21, size=3,aes(fill=Biome))+
  scale_fill_manual(values = col.ecoregion)+
  theme_classic()+xlab('Absolute latitude')+ylab('Linear slope')+
  theme(legend.position = 'none')+
  geom_text_npc(aes(npcx =0.05, npcy =0.75, label='italic(r)==0.482'),parse=TRUE)+
  geom_text_npc(aes(npcx =0.05, npcy =0.65, label='italic(p)<0.001'),parse=TRUE)

# 2. Quadratic slopes vs. Latitude and AI
with(subset(coef.ecoregion, variable=='EcM.p.stem2'), cor.test(Lat_abs, Estimate))
with(subset(coef.ecoregion, variable=='EcM.p.stem2'), cor.test(AI, Estimate))

c.ecoregion=ggplot(subset(coef.ecoregion, variable=='EcM.p.stem2'), aes(Lat_abs, Estimate))+
  geom_smooth(method = glm, col='black',se=F)+
  geom_point(shape=21, size=3,aes(fill=Biome))+
  scale_fill_manual(values = col.ecoregion)+
  theme_classic()+xlab('Absolute latitude')+ylab('Quadratic slope')+
  theme(legend.position = 'none')+
  geom_text_npc(aes(npcx =0.95, npcy =0.15, label='italic(r)==-0.413'),parse=TRUE)+
  geom_text_npc(aes(npcx =0.95, npcy =0.05, label='italic(p)<0.001'),parse=TRUE)



# 3. Linear + Quadratic vs. Latitude and AI
x=subset(coef.ecoregion, variable=='EcM.p.stem') %>%
  mutate(
    linear = Estimate,
    quadratic = subset(coef.ecoregion, variable=='EcM.p.stem2')$Estimate,
    Diff =  linear+quadratic)

x.ecoregion=subset(coef.ecoregion, variable=='EcM.p.stem') %>%
  mutate(
    linear = Estimate,
    quadratic = subset(coef.ecoregion, variable=='EcM.p.stem2')$Estimate,
    Diff =  linear+quadratic) %>%
  mutate(
    hypo='positive unimodal',
    hypo=ifelse(Diff>0&quadratic>=0, 'positive',hypo),
    hypo=ifelse(Diff<0, 'Negtive unimodal',hypo),
    hypo=ifelse(Diff<0&linear<=0, 'Negtive',hypo)
  )

with(x.ecoregion, cor.test(Lat_abs, Diff))
with(x.ecoregion, cor.test(AI, Diff))

diff.ecoregion=ggplot(x.ecoregion, aes(Lat_abs, Diff))+
  geom_smooth(method = glm, col='black',se=F)+
  geom_point(shape=21, size=3,aes(fill=hypo))+
  scale_fill_manual(values = c('#F8766D', '#C77CFF', '#00BFC4', '#7CAE00'))+
  theme_classic()+xlab('Absolute latitude')+ylab('Linear+Quadratic slopes')+
  theme(legend.position = 'none', legend.title = element_blank(), legend.background = element_blank(), plot.title = element_text(face = 'bold'))+ggtitle('b')+
  geom_hline(yintercept = 0, lty=2, col='red')+
  geom_text_npc(aes(npcx =0.95, npcy =0.15, label='italic(r)==0.558'),parse=TRUE)+
  geom_text_npc(aes(npcx =0.95, npcy =0.05, label='italic(p)<0.001'),parse=TRUE)+
  geom_point(data=data.frame(x=rep(22,4), y=seq(8,10, length=4)),aes(x=x, y=y), shape=21,fill=c('#F8766D', '#C77CFF', '#7CAE00', '#00BFC4'),size=2)+
  annotate('text',x=rep(24,4), y=seq(8,10, length=4), label=c('Negative','Negative unimodal','Posive unimodal','Positive'), hjust = 0,size=2)


diff1.ecoregion=ggplot(x.ecoregion, aes(AI, Diff))+
  geom_smooth(method = glm, col='black',se=F)+
  geom_point(shape=21, size=3,aes(fill=hypo))+
  scale_fill_manual(values = c('#F8766D', '#C77CFF', '#00BFC4', '#7CAE00'))+
  theme_classic()+xlab('Aridity index')+ylab('Linear+Quadratic slopes')+
  theme(legend.position = 'none', legend.title = element_blank(), legend.background = element_blank(), plot.title = element_text(face = 'bold'))+ggtitle('c')+
  geom_hline(yintercept = 0, lty=2, col='red')+
  geom_text_npc(aes(npcx =0.95, npcy =0.15, label='italic(r)==-0.165'),parse=TRUE)+
  geom_text_npc(aes(npcx =0.95, npcy =0.05, label='italic(p)==0.114'),parse=TRUE)



# Fig. 3 - Combine plot ---------------------------------------------------
ggsave('figure/Biomes.png', plot=ggarrange(a.biome, biome.legend, diff.biome, diff1.biome, nrow = 2,ncol = 2), width = 6.5, height = 6, dpi = 800)
ggsave('figure/Ecoregion.png', plot=ggarrange(a.ecoregion, ecoregion.legend, diff.ecoregion, diff1.ecoregion,nrow = 2,ncol = 2), width = 6.5, height = 6, dpi = 800)

