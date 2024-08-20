library(piecewiseSEM)
library(tidyverse)

global.pro.sel <- readRDS('data/global.pro.sel.rds')


sem.sr(subset(global.pro.sel), 'EcM.p.stem')
sem.sr(subset(global.pro.sel, region=='Boreal'), 'EcM.p.stem')
sem.sr(subset(global.pro.sel, region=='Temperate'), 'EcM.p.stem')
sem.sr(subset(global.pro.sel, region=='Tropical'), 'EcM.p.stem')


sem.sr(subset(global.pro.sel), 'EcM.p.BA')
sem.sr(subset(global.pro.sel, region=='Boreal'), 'EcM.p.BA')
sem.sr(subset(global.pro.sel, region=='Temperate'), 'EcM.p.BA')
sem.sr(subset(global.pro.sel, region=='Tropical'), 'EcM.p.BA')
