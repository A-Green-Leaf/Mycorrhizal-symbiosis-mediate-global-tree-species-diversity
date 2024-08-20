library(randomForest)
library(tidyverse)
library(ggpubr)
global.pro.sel <- readRDS('data/global.pro.sel.rds')
source('script/functions.R')

global.rf <- rf.j(global.pro.sel)
bore.rf   <- rf.j(subset(global.pro.sel, region=='Boreal'))
temp.rf   <- rf.j(subset(global.pro.sel, region=='Temperate'))
trop.rf   <- rf.j(subset(global.pro.sel, region=='Tropical'))
global.rf; bore.rf; temp.rf; trop.rf

a=plot.imp(global.rf, "Global")
b=plot.imp(bore.rf,   "Boreal")
c=plot.imp(temp.rf,   "Temperate")+ylab('Decrease in node impurities')
d=plot.imp(trop.rf,   "Tropical")+ylab('Decrease in node impurities')
ggsave('figure/Fig.3.png', plot=ggarrange(a,b,c,d, nrow = 2, ncol = 2, labels = c('a','b','c','d')), 
  width=8, height=8,dpi = 600)
