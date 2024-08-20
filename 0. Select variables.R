global.pro.sel <- readRDS('data/global.pro.sel.rds')
global.pro     <- readRDS('data/global.pro.rds')

# First check of all predictors
data.pro.sel <- global.pro.sel %>% 
  select(SR,EcM.p.stem,EcM.p.stem2, bio1,bio12,AI, elev,Slope, sg.BD.0_200:sg.soc.0_200,Area_log, DBH) %>% na.omit()
colnames(data.pro.sel) <- c('SR','EcM','EcM2','MAT','MAP','AI','Elevation','Slope',
  'Soil BD','Soil cec','Soil cfvo','Soil N','Soil pH','Soil clay','Soil sand','Soil silt','Soil soc','Plot size_log','DBH')

png('figure/Fig.S1.cormatrix.png',width = 16,height = 16,units = 'cm',res = 800)
corrplot::corrplot(cor(data.pro.sel),type = "upper", order = "hclust", tl.col = "black", tl.srt = 45)
dev.off()


data.pro <- global.pro %>% 
  select(SR,EcM.p.stem:EcM.p.BA, bio1,bio12, AI:Slope, soilTP, sg.BD:evenness) %>% na.omit()

png('figure/Fig.S1.cormatrix.png',width = 16,height = 16,units = 'cm',res = 800)
corrplot::corrplot(cor(data.pro),type = "upper", order = "hclust", tl.col = "black", tl.srt = 45)
dev.off()

