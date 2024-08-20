library(tidyverse)
library(ggplot2)
library(ggpubr)

# Test hypothesis of species pool change #####
stem.sppool <- readRDS('data/stem.sppool.rds')
abund.sppool <- readRDS('data/abund.sppool.rds')

stem.density <- stem.sppool %>%
  group_by(plot,familyTNRS, genusTNRS, speciesTNRS) %>%
  summarise(
    density_N = length(na.omit(speciesTNRS)),
    Lat = mean(Lat, na.rm=T),
    Lon = mean(Lon, na.rm = T),
    area = mean(area, na.rm = T)
  )

abund.density <- abund.sppool %>%
  mutate(
    density_N = ceiling(density_N)
  )

source('script/functions.R')


col.sel <- c('plot','Lat','Lon','speciesTNRS','genusTNRS','familyTNRS','density_N','area')

all.orginal.data <- rbind(stem.density[,col.sel], abund.density[,col.sel])
all.orginal.data$latbin <- cut(all.orginal.data$Lat, seq(-90, 90, by=10),  labels = seq(-90, 90, by=10)[-1])
all.orginal.data$lonbin <- cut(all.orginal.data$Lon, seq(-180, 180, by=10),labels = seq(-180, 180, by=10)[-1])
all.orginal.data$grid   <- with(all.orginal.data, paste0(latbin, ".", lonbin))
all.orginal.data        <- match.mycor(all.orginal.data)



# Rarefied richness -------------------------------------------------------


all.density <- all.orginal.data %>%
  mutate(mycor=ifelse(mycor=='AM', 'Non_EcM', mycor)) %>%
  group_by(grid, latbin, lonbin, mycor, speciesTNRS) %>%
  summarise(
    density_N = sum(density_N, na.rm = T)
  ) %>% na.omit()

all.density.ecm <- all.density %>% filter(mycor=='EcM')
all.density.nonecm <- all.density %>% filter(mycor=='Non_EcM')

all.density.matrix.ecm <- all.density.ecm %>% 
  as.data.frame() %>%
  select(grid, speciesTNRS, density_N) %>%
  pivot_wider(names_from = speciesTNRS, values_from = density_N) %>%
  mutate(across(-grid, ~replace_na(., 0))) %>%
  column_to_rownames(var = 'grid')

all.density.matrix.nonecm <- all.density.nonecm %>% 
  as.data.frame() %>%
  select(grid, speciesTNRS, density_N) %>%
  pivot_wider(names_from = speciesTNRS, values_from = density_N) %>%
  mutate(across(-grid, ~replace_na(., 0))) %>%
  column_to_rownames(var = 'grid')


data=all.density.matrix.ecm
N=min(rowSums(data))
SR.ecm <- rarefy(data, 100) %>% as.data.frame()
SR.ecm$grid <- rownames(SR.ecm)
colnames(SR.ecm) <- c('SR','grid')
SR.ecm$mycor <- 'EcM'

data=all.density.matrix.nonecm
N=min(rowSums(data))
SR.nonecm <- rarefy(data, 100) %>% as.data.frame()
SR.nonecm$grid <- rownames(SR.nonecm)
colnames(SR.nonecm) <- c('SR','grid')
SR.nonecm$mycor <- 'Non_EcM'

SR.data <- rbind(SR.ecm, SR.nonecm) %>%
  separate(grid, into = c('latbin','lonbin'), sep = '\\.') %>%
  mutate(
    latbin = as.numeric(latbin),
    lonbin = as.numeric(lonbin),
    value = SR
  )


sp.pool <- SR.data
sp.pool$region <- 'Temperate'
sp.pool$region[which(abs(sp.pool$latbin) < 23.5)] <- 'Tropical'
sp.pool$region[which(abs(sp.pool$latbin) >= 50)]  <- 'Boreal'

sp.pool.region <- sp.pool %>% 
  group_by(region, mycor) %>% 
  summarise(mean=mean(value),se=sd(value)/sqrt(length(value)))

a=ggplot(sp.pool, aes(latbin, value, col=mycor,fill=mycor))+
  geom_point(alpha=0.3)+geom_jitter()+theme_classic()+
  xlab('Latitude')+ylab('Richness of species pool')+
  theme(
    legend.position = c(0.8,0.8), 
    legend.title = element_blank(),
    legend.background = element_blank())+
  scale_y_log10()


b=ggplot(sp.pool.region, aes(region, mean,fill=mycor)) + 
  geom_bar(stat = "identity", position = "dodge", alpha = 0.5)  +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), position = position_dodge(0.9), width = 0.25,show.legend = FALSE)+
  theme_classic()+theme(legend.position = 'none')+
  xlab('')+ylab('Richness of species pool')+
  scale_y_log10()



ggsave('figure/Fig.S6.species pool.png', plot=ggarrange(a,b, nrow = 1, ncol = 2,labels = c('a','b')), 
  width=8, height=4,dpi = 600)







