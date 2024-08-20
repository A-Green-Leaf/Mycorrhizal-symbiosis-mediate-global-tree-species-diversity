library(MASS)
library(spdep)
library(tidyverse)
library(ggplot2)
library(ggpubr)

source('script/functions.R')

global.pro.sel <- readRDS('data/global.pro.sel.rds') %>%
  dplyr::select(SR,shannon1,simpson1,
    EcM.p.stem, EcM.p.stem2, EcM.p.BA, EcM.p.BA2,
    bio1,bio12,AI,elev,Slope,sg.clay.0_200,sg.cec.0_200,DBH,Area_log, Lat,Lon,
    Lat_abs, region) %>% na.omit()
global.pro.sel$sg.clay <- global.pro.sel$sg.clay.0_200
global.pro.sel$sg.cec <- global.pro.sel$sg.cec.0_200
# Model 0 - GLM with original predictors without SAC --------------------------------------------

#1. Null model 1
glm.sr.global = glm.nb(SR~(EcM.p.stem+EcM.p.stem2)+
    bio1+AI+elev+Slope+sg.clay+sg.cec+DBH+Area_log, data = global.pro.sel)

#2. Null model 2
glm.sr.globallat = glm.nb(SR~(EcM.p.stem+EcM.p.stem2)+
    Lat_abs+AI+elev+Slope+sg.clay+sg.cec+DBH+Area_log, data = global.pro.sel)

#3. + ecm*lat
glm.sr.lat    = glm.nb(SR~(EcM.p.stem+EcM.p.stem2)*Lat_abs+
    AI+elev+Slope+sg.clay+sg.cec+DBH+Area_log, data = global.pro.sel)

#4. + ecm*ai
glm.sr.ai    = glm.nb(SR~(EcM.p.stem+EcM.p.stem2)*AI+
    Lat_abs+elev+Slope+sg.clay+sg.cec+DBH+Area_log, data = global.pro.sel)

#5. + ecm*lat + ecm*ai
glm.sr.lat.aino    = glm.nb(SR~(EcM.p.stem+EcM.p.stem2)*Lat_abs+(EcM.p.stem+EcM.p.stem2)*AI+
    elev+Slope+sg.clay+sg.cec+DBH+Area_log, data = global.pro.sel)

#6. + ecm*lat*ai
glm.sr.lat.ai    = glm.nb(SR~(EcM.p.stem+EcM.p.stem2)*AI*Lat_abs+
    elev+Slope+sg.clay+sg.cec+DBH+Area_log, data = global.pro.sel)



#1-1. BA model
glm.sr.BA  = glm.nb(SR~(EcM.p.BA+EcM.p.BA2)*AI*Lat_abs+
    elev+Slope+sg.clay+sg.cec+DBH+Area_log, data = global.pro.sel)

#1-2. Shannon model
glm.shannon = glm(shannon1~(EcM.p.stem+EcM.p.stem2)*AI*Lat_abs+
    elev+Slope+sg.clay+sg.cec+DBH+Area_log, data = global.pro.sel, family = Gamma())

#1-3. Simpson model
glm.simpson = glm(simpson1~(EcM.p.stem+EcM.p.stem2)*AI*Lat_abs+
    elev+Slope+sg.clay+sg.cec+DBH+Area_log, data = global.pro.sel, family = Gamma())

# add autocovariate
global.pro.sel$rac.sr.global    <- Spat.cor.rep(glm.sr.global, global.pro.sel[,c('Lat','Lon')], 2000)
global.pro.sel$rac.sr.globallat <- Spat.cor.rep(glm.sr.globallat, global.pro.sel[,c('Lat','Lon')], 2000)
global.pro.sel$rac.sr.lat       <- Spat.cor.rep(glm.sr.lat,    global.pro.sel[,c('Lat','Lon')], 2000)
global.pro.sel$rac.sr.ai        <- Spat.cor.rep(glm.sr.ai,     global.pro.sel[,c('Lat','Lon')], 2000)
global.pro.sel$rac.sr.lat.aino  <- Spat.cor.rep(glm.sr.lat.aino,   global.pro.sel[,c('Lat','Lon')], 2000)
global.pro.sel$rac.sr.lat.ai    <- Spat.cor.rep(glm.sr.lat.ai, global.pro.sel[,c('Lat','Lon')], 2000)

global.pro.sel$rac.sr.BA     <- Spat.cor.rep(glm.sr.BA,     global.pro.sel[,c('Lat','Lon')], 2000)
global.pro.sel$rac.shannon   <- Spat.cor.rep(glm.shannon,   global.pro.sel[,c('Lat','Lon')], 2000)
global.pro.sel$rac.simpson   <- Spat.cor.rep(glm.simpson,   global.pro.sel[,c('Lat','Lon')], 2000)

# Model 1 - GLM with original predictors with SAC --------------------------------------------
                                # Model selection
#1. Null model 1
glm.sr.global1 = glm.nb(SR~(EcM.p.stem+EcM.p.stem2)+
    bio1+AI+elev+Slope+sg.clay+sg.cec+DBH+Area_log+rac.sr.global, data = global.pro.sel)

#2. Null model 2
glm.sr.globallat1 = glm.nb(SR~(EcM.p.stem+EcM.p.stem2)+
    Lat_abs+AI+elev+Slope+sg.clay+sg.cec+DBH+Area_log+rac.sr.globallat, data = global.pro.sel)

#3. + ecm*lat
glm.sr.lat1    = glm.nb(SR~(EcM.p.stem+EcM.p.stem2)*Lat_abs+
    AI+elev+Slope+sg.clay+sg.cec+DBH+Area_log+rac.sr.lat, data = global.pro.sel)

#4. + ecm*ai
glm.sr.ai1    = glm.nb(SR~(EcM.p.stem+EcM.p.stem2)*AI+
    Lat_abs+elev+Slope+sg.clay+sg.cec+DBH+Area_log+rac.sr.ai, data = global.pro.sel)

#5. + ecm*lat + ecm*ai
glm.sr.lat.aino1    = glm.nb(SR~(EcM.p.stem+EcM.p.stem2)*Lat_abs+(EcM.p.stem+EcM.p.stem2)*AI+
    elev+Slope+sg.clay+sg.cec+DBH+Area_log+rac.sr.lat.aino, data = global.pro.sel)

#6. + ecm*lat*ai
glm.sr.lat.ai1    = glm.nb(SR~(EcM.p.stem+EcM.p.stem2)*AI*Lat_abs+
    elev+Slope+sg.clay+sg.cec+DBH+Area_log+rac.sr.lat.ai, data = global.pro.sel)

# glm.sr.lat.ai2 with scaled predictors
# evaluate model performance
data.frame(
  AIC(glm.sr.global1, glm.sr.globallat1, glm.sr.lat1, glm.sr.ai1, glm.sr.lat.aino1, glm.sr.lat.ai1),
  r2=c(get.r2(glm.sr.global1), get.r2(glm.sr.globallat1), get.r2(glm.sr.lat1), get.r2(glm.sr.ai1), get.r2(glm.sr.lat.aino1), get.r2(glm.sr.lat.ai1))
)



#1-1. BA model
glm.sr.BA1  = glm.nb(SR~(EcM.p.BA+EcM.p.BA2)*AI*Lat_abs+
    elev+Slope+sg.clay+sg.cec+DBH+Area_log+rac.sr.BA, data = global.pro.sel)

#1-2. Shannon model
glm.shannon1 = glm(shannon1~(EcM.p.stem+EcM.p.stem2)*AI*Lat_abs+
    elev+Slope+sg.clay+sg.cec+DBH+Area_log+rac.shannon, data = global.pro.sel, family = Gamma())

#1-3. Simpson model
glm.simpson1 = glm(simpson1~(EcM.p.stem+EcM.p.stem2)*AI*Lat_abs+
    elev+Slope+sg.clay+sg.cec+DBH+Area_log+rac.simpson, data = global.pro.sel, family = Gamma())


# Model 2 - GLM with scaled predictors ----------------------------------------------
global.pro.sel.sc <- global.pro.sel
predictor <- c('EcM.p.stem','EcM.p.BA','EcM.p.stem2','EcM.p.BA2','Lat_abs','bio1','AI','elev','Slope', 'sg.clay','sg.cec', 'DBH','Area_log',
  'rac.sr.global','rac.sr.globallat',  'rac.sr.lat','rac.sr.ai',  'rac.sr.lat.aino','rac.sr.lat.ai',  'rac.sr.BA','rac.shannon','rac.simpson')
global.pro.sel.sc <- cbind(global.pro.sel.sc[,c('SR','shannon1','simpson1','region')], scale(global.pro.sel.sc[, predictor]))

#1. Null model 1
glm.sr.global2 = glm.nb(SR~bio1+AI+elev+Slope+sg.clay+sg.cec+DBH+Area_log+rac.sr.global+
    (EcM.p.stem+EcM.p.stem2), data = global.pro.sel.sc)

#2. Null model 2
glm.sr.globallat2 = glm.nb(SR~Lat_abs+AI+elev+Slope+sg.clay+sg.cec+DBH+Area_log+rac.sr.globallat+
    (EcM.p.stem+EcM.p.stem2), data = global.pro.sel.sc)

#3. + ecm*lat
glm.sr.lat2    = glm.nb(SR~AI+elev+Slope+sg.clay+sg.cec+DBH+Area_log+rac.sr.lat+
    (EcM.p.stem+EcM.p.stem2)*Lat_abs, data = global.pro.sel.sc)

#4. + ecm*ai
glm.sr.ai2    = glm.nb(SR~Lat_abs+elev+Slope+sg.clay+sg.cec+DBH+Area_log+rac.sr.ai+
    (EcM.p.stem+EcM.p.stem2)*AI, data = global.pro.sel.sc)

#5. + ecm*lat + ecm*ai
glm.sr.lat.aino2    = glm.nb(SR~elev+Slope+sg.clay+sg.cec+DBH+Area_log+rac.sr.lat.aino+
    (EcM.p.stem+EcM.p.stem2)*Lat_abs+(EcM.p.stem+EcM.p.stem2)*AI, data = global.pro.sel.sc)

#6. + ecm*lat*ai
glm.sr.lat.ai2    = glm.nb(SR~elev+Slope+sg.clay+sg.cec+DBH+Area_log+rac.sr.lat.ai+
    (EcM.p.stem+EcM.p.stem2)*AI*Lat_abs, data = global.pro.sel.sc)


#1-1. BA model
glm.sr.BA2  = glm.nb(SR~elev+Slope+sg.clay+sg.cec+DBH+Area_log+rac.sr.BA+
    (EcM.p.BA+EcM.p.BA2)*AI*Lat_abs, data = global.pro.sel.sc)

#1-2. Shannon model
glm.shannon2 = glm(shannon1~elev+Slope+sg.clay+sg.cec+DBH+Area_log+rac.shannon+
    (EcM.p.stem+EcM.p.stem2)*AI*Lat_abs, data = global.pro.sel.sc, family = Gamma())

#1-3. Simpson model
glm.simpson2 = glm(simpson1~elev+Slope+sg.clay+sg.cec+DBH+Area_log+rac.simpson+
    (EcM.p.stem+EcM.p.stem2)*AI*Lat_abs, data = global.pro.sel.sc, family = Gamma())






# Fig. 2-GLM Lat*AI -----------------------------------------------------
xname <- c('bio1','elev','Slope','sg.clay','sg.cec','DBH','Area_log',
  'rac.sr.global','rac.sr.globallat',  'rac.sr.lat','rac.sr.ai',  'rac.sr.lat.aino','rac.sr.lat.ai',  'rac.sr.BA','rac.shannon','rac.simpson')

pal <- colorRampPalette(c('red','#F0BB41','#619CFF'))

newdata.lat.ai <- expand_grid(
  EcM.p.stem=seq(0,1,length=100),
  Lat_abs = seq(0,60, by=10),
  AI = c(0.5,1.2,2),
  as.data.frame(t(apply(global.pro.sel[, xname], 2, median)))
) %>% 
  mutate(
    EcM.p.stem2 = EcM.p.stem^2,
    EcM.p.BA = EcM.p.stem,
    EcM.p.BA2 = EcM.p.stem2)

newdata.lat.ai$pred         <- predict(glm.sr.lat.ai1, newdata.lat.ai, type = "response", se.fit=TRUE)$fit
newdata.lat.ai$pred.ba      <- predict(glm.sr.BA1,     newdata.lat.ai, type = "response", se.fit=TRUE)$fit
newdata.lat.ai$pred.shannon <- predict(glm.shannon1,   newdata.lat.ai, type = "response", se.fit=TRUE)$fit
newdata.lat.ai$pred.simpson <- predict(glm.simpson1,   newdata.lat.ai, type = "response", se.fit=TRUE)$fit

p.forest=forestplot(glm.sr.lat.ai2, ' ')+
  theme(axis.text.x = element_text(angle = 15, hjust=1),
    plot.margin = ggplot2::margin(5,5,5,10))
p1 <- ggarrange(p.forest, nrow = 1, ncol = 1, labels = 'a')

p <- ggplot()+
  theme_classic()+
  scale_y_continuous(trans='log10')+
  ggtitle(' ')+
  xlab('EcM stem proportion')+
  ylab('Species richness')+
  theme(plot.title = element_text(hjust = 0.5),
    legend.background = element_blank(),
    legend.position = c(0.83,0.75),
    legend.key.height = unit(0.15,'cm'))

a=p+geom_line(data = subset(newdata.lat.ai, AI==0.5), aes(EcM.p.stem, pred, color=factor(Lat_abs)), linewidth=0.8)+
  scale_color_manual(values = pal(7))+
  labs(color='Latitude ( )')+ggtitle('High aridity')
b=p+geom_line(data = subset(newdata.lat.ai, AI==1.2), aes(EcM.p.stem, pred, color=factor(Lat_abs)), linewidth=0.8)+
  scale_color_manual(values = pal(7))+
  theme(legend.position = 'none')+ggtitle('Middle aridity')
c=p+geom_line(data = subset(newdata.lat.ai, AI==2), aes(EcM.p.stem, pred, color=factor(Lat_abs)), linewidth=0.8)+
  scale_color_manual(values = pal(7))+
  theme(legend.position = 'none')+ggtitle('Low aridity')

p2 <- ggarrange(a,b,c, nrow = 1, ncol = 3, labels = c('b','c','d'))

ggsave('figure/model5.png', plot = ggarrange(p1, p2, nrow = 2,ncol = 1), 
  height = 5.5, width = 9, dpi = 800)



# Fig. S4-GLM BA, Shannon, Simpson ----------------------------------------
p <- ggplot()+
  theme_classic()+
  geom_point(data = global.pro.sel, aes(EcM.p.stem, SR), col='gray',shape=20,size=2,alpha=0.3)+
  scale_y_continuous(trans='log10')+
  ggtitle(' ')+
  xlab('EcM BA proportion')+
  ylab('Species richness')+
  theme(plot.title = element_text(face = 'bold'),
    legend.background = element_blank(),
    legend.position = c(0.75,0.85),
    legend.key.height = unit(0.3,'cm'))

a1=p+geom_line(data = subset(newdata.lat.ai, AI==0.5), aes(EcM.p.stem, pred.ba, color=factor(Lat_abs)), linewidth=0.8)+
  scale_color_manual(values = pal(7))+
  labs(color='Absolute latitude ( )')+ggtitle('a')
a2=p+geom_line(data = subset(newdata.lat.ai, AI==1.2), aes(EcM.p.stem, pred.ba, color=factor(Lat_abs)), linewidth=0.8)+
  scale_color_manual(values = pal(7))+
  theme(legend.position = 'none')+ggtitle(' ')
a3=p+geom_line(data = subset(newdata.lat.ai, AI==2.0), aes(EcM.p.stem, pred.ba, color=factor(Lat_abs)), linewidth=0.8)+
  scale_color_manual(values = pal(7))+
  theme(legend.position = 'none')+ggtitle(' ')

p <- ggplot()+
  theme_classic()+
  geom_point(data = global.pro.sel, aes(EcM.p.stem, shannon1), col='gray',shape=20,size=2,alpha=0.3)+
  ggtitle(' ')+
  xlab('EcM stem proportion')+
  ylab('Shannon')+
  theme(plot.title = element_text(face = 'bold'),
    legend.background = element_blank(),
    legend.position = c(0.75,0.85),
    legend.key.height = unit(0.3,'cm'))

b1=p+geom_line(data = subset(newdata.lat.ai, AI==0.5), aes(EcM.p.stem, pred.shannon, color=factor(Lat_abs)), linewidth=0.8)+
  scale_color_manual(values = pal(7))+
  labs(color='Absolute latitude ( )')+ggtitle('b')
b2=p+geom_line(data = subset(newdata.lat.ai, AI==1.2), aes(EcM.p.stem, pred.shannon, color=factor(Lat_abs)), linewidth=0.8)+
  scale_color_manual(values = pal(7))+
  theme(legend.position = 'none')+ggtitle(' ')
b3=p+geom_line(data = subset(newdata.lat.ai, AI==2.0), aes(EcM.p.stem, pred.shannon, color=factor(Lat_abs)), linewidth=0.8)+
  scale_color_manual(values = pal(7))+
  theme(legend.position = 'none')+ggtitle(' ')


p <- ggplot()+
  theme_classic()+
  geom_point(data = global.pro.sel, aes(EcM.p.stem, simpson1), col='gray',shape=20,size=2,alpha=0.3)+
  ggtitle(' ')+
  xlab('EcM stem proportion')+
  ylab('Simpson')+
  theme(plot.title = element_text(face = 'bold'),
    legend.background = element_blank(),
    legend.position = c(0.75,0.85),
    legend.key.height = unit(0.3,'cm'))

c1=p+geom_line(data = subset(newdata.lat.ai, AI==0.5), aes(EcM.p.stem, pred.simpson, color=factor(Lat_abs)), linewidth=0.8)+
  scale_color_manual(values = pal(7))+
  labs(color='Absolute latitude ( )')+ggtitle('c')
c2=p+geom_line(data = subset(newdata.lat.ai, AI==1.2), aes(EcM.p.stem, pred.simpson, color=factor(Lat_abs)), linewidth=0.8)+
  scale_color_manual(values = pal(7))+
  theme(legend.position = 'none')+ggtitle(' ')
c3=p+geom_line(data = subset(newdata.lat.ai, AI==2.0), aes(EcM.p.stem, pred.simpson, color=factor(Lat_abs)), linewidth=0.8)+
  scale_color_manual(values = pal(7))+
  theme(legend.position = 'none')+ggtitle(' ')


ggsave('figure/Fig.S4.png', 
  plot = ggarrange(a1,a2,a3,b1,b2,b3,c1,c2,c3, nrow = 3, ncol = 3), height = 9, width = 9)






# Fig. S3 by three regions ####

# Boreal
data.region <- subset(global.pro.sel, region=='Boreal')
mod <- MASS::glm.nb(SR~EcM.p.stem+EcM.p.stem2, data = data.region)
new.data <- data.frame(EcM.p.stem=seq(0,1, length=100)) %>% mutate(EcM.p.stem2 = EcM.p.stem^2)
new.data$pred <- predict(mod, new.data, type = "response", se.fit=TRUE)$fit
new.data$se <- predict(mod, new.data, type = "response", se.fit=TRUE)$se.fit

a=ggplot()+theme_classic()+
  geom_point(data=data.region,aes(EcM.p.stem, SR),col="gray",shape=20,size=2,alpha=0.3)+
  xlab("EcM stem proportion")+ylab("Species richness")+ggtitle('Boreal')+
  geom_ribbon(data=new.data, aes(x=EcM.p.stem, y=pred, ymin = pred-se, ymax = pred+se, fill = 'cyan4'), alpha = .5)+
  geom_line(data=new.data, aes(EcM.p.stem, pred), col='cyan4', size = 1)+
  scale_y_continuous(trans='log10')+
  theme(legend.position = 'none')


# Temperate
data.region <- subset(global.pro.sel, region=='Temperate')
mod <- MASS::glm.nb(SR~EcM.p.stem+EcM.p.stem2, data = data.region)
new.data <- data.frame(EcM.p.stem=seq(0,1, length=100)) %>% mutate(EcM.p.stem2 = EcM.p.stem^2)
new.data$pred <- predict(mod, new.data, type = "response", se.fit=TRUE)$fit
new.data$se <- predict(mod, new.data, type = "response", se.fit=TRUE)$se.fit

b=ggplot()+theme_classic()+
  geom_point(data=data.region,aes(EcM.p.stem, SR),col="gray",shape=20,size=2,alpha=0.3)+
  xlab("EcM stem proportion")+ylab("Species richness")+ggtitle('Temperate')+
  geom_ribbon(data=new.data, aes(x=EcM.p.stem, y=pred, ymin = pred-se, ymax = pred+se, fill = 'cyan4'), alpha = .15)+
  geom_line(data=new.data, aes(EcM.p.stem, pred), col='cyan4', size = 1)+
  scale_y_continuous(trans='log10')+
  theme(legend.position = 'none')


# Tropical
data.region <- subset(global.pro.sel, region=='Tropical')
mod <- MASS::glm.nb(SR~EcM.p.stem+EcM.p.stem2, data = data.region)
new.data <- data.frame(EcM.p.stem=seq(0,1, length=100)) %>% mutate(EcM.p.stem2 = EcM.p.stem^2)
new.data$pred <- predict(mod, new.data, type = "response", se.fit=TRUE)$fit
new.data$se <- predict(mod, new.data, type = "response", se.fit=TRUE)$se.fit

c=ggplot()+theme_classic()+
  geom_point(data=data.region,aes(EcM.p.stem, SR),col="gray",shape=20,size=2,alpha=0.3)+
  xlab("EcM stem proportion")+ylab("Species richness")+ggtitle('Tropical')+
  geom_ribbon(data=new.data, aes(x=EcM.p.stem, y=pred, ymin = pred-se, ymax = pred+se, fill = 'cyan4'), alpha = .15)+
  geom_line(data=new.data, aes(EcM.p.stem, pred), col='cyan4', size = 1)+
  scale_y_continuous(trans='log10')+
  theme(legend.position = 'none')

ggsave('figure/Fig.S3.png', plot = ggarrange(a,b,c,nrow = 1,ncol = 3, labels = c('a','b','c')), width = 10, height = 3.3, dpi = 800)




# Table S2--Coefficients ####
glm.sr.null <- glm.nb(SR~EcM.p.stem+EcM.p.stem2, data = global.pro.sel)

table.s2.coef <- rbind(
  glm.coef.sc(glm.sr.null, 'Null'), 
  glm.coef.sc(glm.sr.lat.ai2, 'Latitude * AI'), 
  glm.coef.sc(glm.sr.BA2,   'Basal area'), 
  glm.coef.sc(glm.shannon2,  'Shannon'), 
  glm.coef.sc(glm.simpson2, 'Simpson')) %>%
  select(mod.name, variable, Estimate, SE, p, lower, upper) %>%
  mutate(
    Estimate = round(Estimate, 3),
    SE = round(SE, 3),
    p = round(p, 4),
    lower = round(lower, 3),
    upper = round(upper, 3))

write.csv(table.s2.coef, 'table/Table.S2.coef.csv')

# Table S3--VIF ####
table.s3.vif <- rbind( 
  get.glm.vif(glm.sr.lat.ai2, 'Latitude * AI'), 
  get.glm.vif(glm.sr.BA2, 'Basal area'), 
  get.glm.vif(glm.shannon2,  'Shannon'), 
  get.glm.vif(glm.simpson2, 'Simpson')) 

table.s3.vif <- subset(table.s3.vif, Variable!='Lat_abs'&Variable!='AI')

write.csv(table.s3.vif, 'table/Table.s3.vif.csv')


