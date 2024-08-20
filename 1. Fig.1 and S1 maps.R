library(ggplot2)
library(sf)
library(ggpubr)

global.pro     <- readRDS('data/global.pro.rds')
global.pro.sel <- readRDS('data/global.pro.sel.rds')
ecoregion <- st_read("D:/0. Publised dataset/2. Environments/5. Ecoregion/Olson-Dinerstein-2017-Ecoregions/Ecoregions2017.shp")[,c('ECO_NAME','BIOME_NAME','REALM')]


# Fig.1-Hypothesis and global map -------------------------------------------------------------------

# Set biome color
levels(factor(unique(ecoregion$BIOME_NAME)))
bor='darkseagreen1';tem='darkolivegreen3';tro='darkgreen';ot='lightgray'
fill.cor = c(bor, ot,ot,ot,ot,ot,ot, tem,tem, ot, tro,tro, ot, tro, ot)

y=alpha('darkolivegreen3', 0.5);ot='lightgray'
fill.cor = c(y, ot,ot,ot,ot,ot,ot, y,y, ot, y,y, ot, y, ot)



# Fig.1a-1c
x=seq(0,1,length=100)
#y1=-100*x+101; y2=100*x-100*x^2+1; y3=100*x-120*x^2+21
y1= (-1*x+1)*100; y2= (1*x-1*x^2)*400; y3=(1*x-1.5*x^2+0.5)*150
p1=ggplot()+
  geom_line(aes(x=x, y=y1),size=1.5)+
  annotate("text",0.5,15,label="y == f(-x)",parse=T,size=6)+
  theme_bw()+theme(plot.title = element_text(hjust = 0.5))+
  labs(title = 'EcM dominance hypothesis', x='', y='Species richness')
p2=ggplot()+
  geom_line(aes(x=x, y=y2),size=1.5)+
  annotate("text",0.5,15,label="y == f(x - x^2)",parse=T,size=6)+
  theme_bw()+theme(plot.title = element_text(hjust = 0.5))+
  labs(title = 'Mycorrhizal mixture hypothesis', x='EcM proportion', y='')
p3=ggplot()+
  geom_line(aes(x=x, y=y3),size=1.5)+
  annotate("text",0.5,15,label="y == f(x - 1.5*x^2)",parse=T,size=6)+
  theme_bw()+theme(plot.title = element_text(hjust = 0.5))+
  labs(title = 'Integrated hypothesis')+xlab('')+ylab('')
p= ggarrange(p1,p2,p3, nrow = 1, ncol = 3, labels = c('a','b','c'))


# Fig.1d
map.sel=ggplot() + 
  geom_map(data = map_data("world"), map = map_data("world"), aes(long, lat, map_id = region), fill=alpha("lightgrey",0.8)) +
  geom_sf(data = ecoregion, aes(fill=BIOME_NAME),color=NA, show.legend = FALSE)+
  scale_fill_manual(values = fill.cor)+
  theme_classic()+
  theme(
    legend.position = 'bottom',
    legend.key.width = unit(2, 'cm'),
    legend.key.height = unit(0.2, 'cm'),
    legend.box.margin = ggplot2::margin(-30,-0,-0,-10),
    legend.spacing.x = unit(0.5, 'cm'),
    axis.line.x = element_blank(),
    axis.line.y = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks = element_blank(),
    plot.margin = ggplot2::margin(0,2,2,2))+
  geom_point(data=global.pro.sel[,c('Lon','Lat','EcM.p.stem')], aes(Lon, Lat, col=EcM.p.stem), size=0.2) +
  scale_color_gradient2(low ='blue',high ='red',midpoint =0.5)+
  labs(color="EcM stem proportion", x="", y="") 

# Fig.1e
mod <- MASS::glm.nb(SR~EcM.p.stem+EcM.p.stem2, data = global.pro.sel)
new.data <- data.frame(EcM.p.stem=seq(0,1, length=100)) %>%
  mutate(EcM.p.stem2 = EcM.p.stem^2)
new.data$pred <- predict(mod, new.data, type = "response", se.fit=TRUE)$fit
new.data$se <- predict(mod, new.data, type = "response", se.fit=TRUE)$se.fit

b<-ggplot()+#ylim(0,300)+
  geom_point(data=global.pro.sel,aes(EcM.p.stem, SR),col="gray",shape=20,size=2,alpha=0.3)+
  xlab("EcM stem proportion")+ylab("Species richness")+
  geom_ribbon(data=new.data, aes(x=EcM.p.stem, y=pred, ymin = pred-se, ymax = pred+se, fill = 'cyan4'), alpha = .15)+
  geom_line(data=new.data, aes(EcM.p.stem, pred), col='cyan4', size = 1)+
  scale_y_continuous(trans='log10')+
  theme_classic()+
  theme(
    plot.margin = ggplot2::margin(7,4,4,4),
    legend.position = 'none')
  
p1 <- ggarrange(map.sel, b, nrow = 1, ncol = 2, widths = c(0.6,0.4), labels = c('d','e'))

ggsave('figure/Fig.1v2.png', plot=ggarrange(p, p1, nrow = 2, ncol = 1, heights = c(0.45, 0.55)),
  width=10.5, height=6,dpi=800)




# Fig. S1-Map for all plots -----------------------------------------------------------------

# Fig.S1 map with all plots
a=ggplot() + 
  geom_map(data = map_data("world"), map = map_data("world"), aes(long, lat, map_id = region), fill=alpha("lightgrey",0.8)) +
  geom_sf(data = ecoregion, aes(fill=BIOME_NAME),color=NA, show.legend = FALSE)+
  scale_fill_manual(values = fill.cor)+
  theme_classic()+
  theme(
    legend.position = 'bottom',
    legend.key.width = unit(3.5, 'cm'),
    legend.key.height = unit(0.5, 'cm'),
    legend.box.margin = ggplot2::margin(-30,-30,-0,-30),
    axis.line.x = element_blank(),
    axis.line.y = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks = element_blank(),
    plot.margin = ggplot2::margin(-1,-1,-1,1),
    plot.title = element_text(vjust = -4,hjust = 0.1))+
  geom_point(data=global.pro[,c('Lon','Lat','EcM.p.stem')], aes(Lon, Lat, col=EcM.p.stem), size=0.3) +
  scale_color_gradient2(low ='blue',high ='red',midpoint =0.5)+
  labs(color="EcM stem proportion",title='Distribution of all forest plots (N=442384)', x="", y="") 

ggsave('figure/Fig.S1-map with all plots.png', plot=a, width=10, height=5,dpi=800)



# Fig. S - AM and diversity -------------------------------------------------
a=ggplot(global.pro.sel, aes(EcM.p.stem, AM.p.stem))+
  geom_point(col="gray",shape=20,size=2,alpha=0.6)+
  xlab('EcM stem proportion')+ylab('AM stem proportion')+
  theme_classic()

# AM vs. richness
# Fig.1e
mod <- MASS::glm.nb(SR~AM.p.stem+AM.p.stem2, data = global.pro.sel)
new.data <- data.frame(AM.p.stem=seq(0,1, length=100)) %>%
  mutate(AM.p.stem2 = AM.p.stem^2)
new.data$pred <- predict(mod, new.data, type = "response", se.fit=TRUE)$fit
new.data$se <- predict(mod, new.data, type = "response", se.fit=TRUE)$se.fit

b<-ggplot()+#ylim(0,300)+
  geom_point(data=global.pro.sel,aes(AM.p.stem, SR),col="gray",shape=20,size=2,alpha=0.3)+
  xlab("AM stem proportion")+ylab("Species richness")+
  geom_ribbon(data=new.data, aes(x=AM.p.stem, y=pred, ymin = pred-se, ymax = pred+se, fill = 'cyan4'), alpha = .15)+
  geom_line(data=new.data, aes(AM.p.stem, pred), col='cyan4', size = 1)+
  scale_y_continuous(trans='log10')+
  theme_classic()+
  theme(
    plot.margin = ggplot2::margin(7,4,4,4),
    legend.position = 'none')

ggsave('figure/Fig.AM and richness.png', plot = ggarrange(a,b,nrow = 1,ncol = 2,labels = c('a','b')),
  width = 6, height = 3, dpi = 800)



