# Assign mycor data function ..... must have species, genus, and family names
match.mycor  <- function(data=china, sp="speciesTNRS", genus="genusTNRS", family='familyTNRS'){
  
  mycor.sp    <- openxlsx::read.xlsx("D:/0. Publised dataset/3. Functional traits/6. FungalRoot/Mycorrhizaltype_data.xlsx", sheet = "Final-MT-Species")
  mycor.genus <- read.csv("D:/0. Publised dataset/3. Functional traits/6. FungalRoot/FungalRoot.csv"); 
  mycor.genus$mycor.genus <- mycor.genus$Mycor; 
  mycor.genus <- mycor.genus[,-2]
  
  mycor.family <- data.frame(
    Family = c('Pinaceae','Betulaceae','Fagaceae','Dipterocarpaceae'),
    mycor.family = 'EcM'
  )
  
  data = merge(data, mycor.sp,    by.x=sp,     by.y="Species",all.x=T,sort=F) # match using species-level data
  data = merge(data, mycor.genus, by.x=genus,  by.y="Genus",  all.x=T,sort=F) # match using genus-level data
  data = merge(data, mycor.family,by.x=family, by.y="Family", all.x=T,sort=F) # match using family-level data
  
  # set mycor using species, genus, family level mycor
  data$mycor = with(data, ifelse(is.na(mycor.genus), Mycor.sp,     mycor.genus)) # genus first, species then
  data$mycor = with(data, ifelse(is.na(mycor),       mycor.family, mycor)) # match family records
  
  # for EcM-AM genus, we should assign to species level
  data$mycor = with(data, ifelse(mycor.genus=='EcM-AM' & !is.na(Mycor.sp) & mycor.genus != Mycor.sp,       Mycor.sp, mycor))
  
  data$mycor[is.na(data$mycor)] = 'unknown' # unknown
  data$mycor[data$mycor=="uncertain"] <- "unknown"
  data$mycor <- gsub("-", "_", data$mycor)
  
  
  data$mycor.origin <- data$mycor # store the origin definition
  
  # Code to replace other mycor to Non-EcM
  data$mycor[-which(data$mycor %in% c('EcM','EcM_AM','AM'))] <- "Non_EcM" # define Non_EcM species
  
  data = subset(data, select = -c(Mycor.sp, mycor.genus, mycor.family))
  return(data)
}

# RAC function - autovariate
Spat.cor <- function(mod, dat, dist) {
  coords <- cbind(dat$Lon, dat$Lat)
  matrix.dist = as.matrix(dist(cbind(dat$Lon, dat$Lat)))
  matrix.dist[1:10, 1:10]
  matrix.dist.inv <- 1/matrix.dist
  matrix.dist.inv[1:10, 1:10]
  diag(matrix.dist.inv) <- 0
  matrix.dist.inv[1:10, 1:10]
  myDist = dist
  rac <- autocov_dist(resid(mod), coords, nbs = myDist, type = "inverse", zero.policy = TRUE, style = "W", longlat = T)
  return(rac)
}

#RAC function when locations repeat (shift latlon)
Spat.cor.rep <- function(mod, dat, dist) {
  coords <- cbind(dat$Lon, dat$Lat) + matrix(runif(2*nrow(dat), 0, 0.00001), nrow = nrow(dat), ncol = 2)
  matrix.dist = as.matrix(dist(cbind(dat$Lon, dat$Lat)))
  matrix.dist[1:10, 1:10]
  matrix.dist.inv <- 1/matrix.dist
  matrix.dist.inv[1:10, 1:10]
  diag(matrix.dist.inv) <- 0
  matrix.dist.inv[1:10, 1:10]
  myDist = dist
  rac <- autocov_dist(resid(mod), coords, nbs = myDist, type = "inverse", zero.policy = TRUE, style = "W", longlat = T)
  return(rac)
}

# Function to forest plot with scaled coefficients
glm.coef.sc <- function(mod=glm.sr.lat2, mod.name='Global'){
  
  output <- summary(mod)$coef[-1,] %>% as.data.frame()
  colnames(output) <- c('Estimate','SE','stat.value','p')
  output$var.origin <- rownames(output)
  #output$p.aov <- anova(mod)$`Pr(>Chi)`[-1]
  
  output$sh <- ifelse(output$p<0.05, alpha('gray',0.7), 'white')
  #output$sh.aov <- ifelse(output$p.aov<0.05, alpha('gray',0.7), 'white')
  
  # confidence interval
  conf <- confint(mod)[-1,]
  colnames(conf) <- c('lower', 'upper')
  
  output <- cbind(output, conf)
  output$mod.name <- mod.name
  output$sample.size <- nrow(mod$model)
  
  # rename variables
  output$variable <- output$var.origin
  output$variable <- str_replace_all(output$variable, '.*rac.*', 'Spatial autocorrelation')
  output$variable <- gsub('EcM.p.stem', 'EcM', output$variable)
  output$variable <- gsub('EcM.p.BA', 'EcM', output$variable)
  
  rename_vector <- c("elev" = "Elevation", "sg.clay" = "Soil clay", 'sg.cec'='Soil cec')
  output <- output %>% 
    mutate(variable=dplyr::recode(variable, !!!rename_vector),
    variable = factor(variable, levels=variable))
  
  return(output)
}

# obtain vif values for glm models
get.glm.vif <- function(mod=glm.sr.region2, mod.name='Global'){
  vif.data <- as.data.frame(car::vif(mod))
  
  if(ncol(vif.data)==3){
    vif.data <- vif.data[,-c(1,2),drop=F]  
  }
  
  vif.data$Variable <- rownames(vif.data)
  vif.data$Model <- mod.name
  
  colnames(vif.data)[1] <- 'vif'
  vif.data$vif <- round(vif.data$vif, 1)
  
  vif.data <- vif.data %>%
    filter(Variable %in% c('elev','Slope','sg.clay','sg.cec','DBH','Area_log','region',
      'rac.sr.global','rac.sr.lat','rac.sr.lat.ai','rac.sr.region','rac.sr.BA','rac.shannon','rac.simpson'))
  
  return(vif.data[,c('Model','Variable','vif')])
}

# Function to plot error plot
forestplot <- function(mod=glm.sr.region2, title='Global'){
  
  output <- glm.coef.sc(mod) %>%
    filter(variable != 'Spatial autocorrelation')

  N1= which(output$variable =='EcM')-1
  N2= nrow(output)
  shade <- data.frame(
    xmin = c(0.7,N1+0.7), xmax = c(N1+0.3,N2+0.3),  
    ymin = - Inf,   ymax = Inf)
  
  ggplot(data=output) +
    geom_point(aes(x=variable, y=Estimate), col="white")+
    geom_rect(data=shade, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill=c('gray','coral'), alpha=0.2)+
    geom_errorbar(aes(x = variable, ymin = lower, ymax = upper), width = 0, linewidth = 0.8)+
    geom_point(aes(x=variable, y=Estimate),fill=output$sh, size = 3, stroke = 1, shape = 21)+
     theme_classic() + ylab("Coefficient estimate") + xlab("") + #coord_flip() +
    geom_hline(yintercept = 0,lty=2) +
    ggtitle(title) +  
    theme(legend.position = "none", plot.title = element_text(face = 'bold'))
  
}

get.r2 <- function(mod){with(mod, 1-(deviance/null.deviance))}


glm.biome <- function(data, name){
  # mod
  mod = glm(SR~EcM.p.stem+EcM.p.stem2+
      bio1+AI+elev+Slope+sg.clay.0_200+sg.cec.0_200+DBH+Area_log, data = data, family = poisson())
  
  # coefficients
  coef <- as.data.frame(summary(mod)$coefficients)
  coef$name <- name
  coef$Lat_abs <- median(data$Lat_abs)
  coef$AI <- median(data$AI)
  coef$variable <- rownames(coef)
  coef$vif <- c(NA, car::vif(mod))
  coef$SR <- median(data$SR)
  
  # r2 and overdispersion
  r2 <- with(mod, (null.deviance-deviance)/null.deviance)
  overdis <- round(with(mod, deviance/df.residual),2)
  
  if(overdis >1) {
    mod = glm.nb(SR~EcM.p.stem+EcM.p.stem2+
        bio1+AI+elev+Slope+sg.clay.0_200+sg.cec.0_200+DBH+Area_log, data = data)
    
    overdis <- round(with(mod, deviance/df.residual),2)
  }
  
  newdata <- expand_grid(
    EcM.p.stem=seq(0,1,length=100),
    as.data.frame(t(apply(data[, 4:11], 2, median)))) %>%
    mutate(EcM.p.stem2 = EcM.p.stem^2)
  
  newdata$pred <- predict(mod, newdata, type = "response", se.fit=TRUE)$fit
  newdata$name <- name
  newdata <- subset(newdata, pred>=1)
  output <- list(mod=mod, coef=coef, r2=r2, overdis=overdis, newdata=newdata)
  
  return(output)
}

glm.ecoregion <- function(data, name){
  # mod
  mod = glm(SR~EcM.p.stem+EcM.p.stem2+
      bio1+AI+elev+Slope+sg.clay.0_200+sg.cec.0_200+DBH+Area_log, data = data, family = poisson())
  
  # coefficients
  coef <- as.data.frame(summary(mod)$coefficients)
  coef$name <- name
  coef$Lat_abs <- median(data$Lat_abs)
  coef$AI <- median(data$AI)
  coef$variable <- rownames(coef)
  coef$vif <- c(NA, vif(mod))
  coef$SR <- median(data$SR)
  
  # r2 and overdispersion
  r2 <- with(mod, (null.deviance-deviance)/null.deviance)
  overdis <- round(with(mod, deviance/df.residual),2)
  
  if(overdis >1) {
    mod = glm.nb(SR~EcM.p.stem+EcM.p.stem2+
        bio1+AI+elev+Slope+sg.clay.0_200+sg.cec.0_200+DBH+Area_log, data = data)
    
    overdis <- round(with(mod, deviance/df.residual),2)
  }
  
  # predict
  predictor <- data %>% dplyr::select(bio1:Area_log)
  newdata <- expand_grid(
    EcM.p.stem=seq(0,1,length=100),
    as.data.frame(t(apply(predictor, 2, median)))) %>%
    mutate(EcM.p.stem2 = EcM.p.stem^2)
  
  newdata$pred <- predict(mod, newdata, type = "response", se.fit=TRUE)$fit
  newdata$name <- name
  newdata <- subset(newdata, pred>=1)
  output <- list(mod=mod, coef=coef, r2=r2, overdis=overdis, newdata=newdata)
  
  return(output)
}



rf.j <- function(data){
  
  all.var <- c('SR',  
    'EcM.p.stem', 
    'bio1','AI',
    'elev','Slope',  
    'sg.clay.0_200', 'sg.cec.0_200',
    'DBH','Area_log')
  data <- na.omit(data[,all.var])
  colnames(data) <- c('SR',  
    'EcM_proportion', 
    'MAT','AI',
    'Elevation','Slope',  
    'Soil_clay', 'Soil_cec',
    'DBH','Plot_size_log')
  
  y=data$SR
  x=data[,-1]
  
  output <- randomForest(x=x,y=y, importance=T)
  return(output)
  
}

# function for IncNodePurity with rank
imp.extr <- function(M){
  res <- data.frame(var = rownames(importance(M)), # variable name
    importance(M),
    imp=importance(M)[,2])
  
  res      <- res[order(res$imp, decreasing = T),]
  res$var1 <- factor(res$var, levels = res$var)
  res$col  <- 'a'
  res$col[which(res$var1 %in% c('EcM_proportion',  'EcM_proportion2'))]   <- 'b'
  #res$col[which(res$var1 %in% c('EcM_stem','EcM_stem2'))] <- 'c'
  #res$col[which(res$var1 %in% c('EcM_BA',  'EcM_BA2'))]   <- 'd'
  
  return(res)
} 

# function for plot relative importance
plot.imp <- function(mod, title){
  ggplot(imp.extr(mod), aes(var1, imp,fill=col)) +  
    geom_col(alpha=0.5,col=alpha('black',0.7)) +  
    theme_classic() + xlab("") + ylab("") + coord_flip() +
    scale_fill_manual(values=c("grey70",'coral')) +
    ggtitle(title)+
    theme(legend.position="none", plot.title = element_text(hjust = 0.5))
}



# Funct. EcM prop. function for stem data
get.pro.stem <- function(data){
  library(vegan)
  
  data$BA   <- with(data, pi*(dbh.cm/2)*(dbh.cm/2))
  # Mycor
  var.name  <- c('plot','speciesTNRS', 'genusTNRS', 'familyTNRS','BA')
  data      <- data[, var.name]
  data$genus<- gsub( " .*$", "", data$speciesTNRS)
  data      <- match.mycor(data, 'speciesTNRS', 'genusTNRS', 'familyTNRS')
  
  # Diversity
  data.abund <- xtabs(data = data,  ~plot+speciesTNRS)  # abundance matrix
  data.div <- data.frame(
    plot    = rownames(data.abund),  
    SR      = specnumber(data.abund), 
    shannon = diversity(data.abund, "shannon"), 
    simpson = diversity(data.abund, "simpson"))
  
  # EcM pro
  abund.prop  <- data %>% 
    group_by(plot, mycor) %>% 
    summarise(
      SR   = length(unique(na.omit(speciesTNRS))), 
      stem = length(na.omit(speciesTNRS)),
      BA   = sum(BA, na.rm=T))
  
  abund.prop  <- maditr::dcast(data.table::setDT(abund.prop), plot~mycor, value.var = c("SR", "stem","BA"))
  
  abund.prop$total.SR   <- rowSums(abund.prop[, paste0("SR_",c('AM', "EcM_AM", "EcM", "Non_EcM"))], na.rm = T)
  abund.prop$EAM.p.SR   <- with(abund.prop, rowSums(cbind(SR_EcM_AM), na.rm = T)/total.SR)
  abund.prop$EcM.p.SR   <- with(abund.prop, rowSums(cbind(SR_EcM, SR_EcM_AM/2),na.rm = T)/total.SR)
  abund.prop$AM.p.SR    <- with(abund.prop, rowSums(cbind(SR_AM, SR_EcM_AM/2),na.rm = T)/total.SR)
  
  abund.prop$total.stem <- rowSums(abund.prop[, paste0("stem_",c('AM', "EcM_AM", "EcM", "Non_EcM"))], na.rm = T)
  abund.prop$EAM.p.stem <- with(abund.prop, rowSums(cbind(stem_EcM_AM), na.rm = T)/total.stem)
  abund.prop$EcM.p.stem <- with(abund.prop, rowSums(cbind(stem_EcM, stem_EcM_AM/2),na.rm = T)/total.stem)
  abund.prop$AM.p.stem  <- with(abund.prop, rowSums(cbind(stem_AM, stem_EcM_AM/2),na.rm = T)/total.stem)
  
  abund.prop$total.BA   <- rowSums(abund.prop[, paste0("BA_",c('AM', "EcM_AM", "EcM", "Non_EcM"))], na.rm = T)
  abund.prop$EAM.p.BA   <- with(abund.prop, rowSums(cbind(BA_EcM_AM), na.rm = T)/total.BA)
  abund.prop$EcM.p.BA   <- with(abund.prop, rowSums(cbind(BA_EcM, BA_EcM_AM/2),na.rm = T)/total.BA)
  abund.prop$AM.p.BA    <- with(abund.prop, rowSums(cbind(BA_AM, BA_EcM_AM/2),na.rm = T)/total.BA)
  
  abund.prop <- merge(abund.prop, data.div, by="plot", all.x=T)
  
  abund.prop <- abund.prop[,c('plot','total.SR','SR','shannon','simpson',
    'EcM.p.SR','EcM.p.stem','EcM.p.BA','AM.p.SR','AM.p.stem','AM.p.BA','EAM.p.SR','EAM.p.stem','EAM.p.BA')]
  
  abund.prop <- as.data.frame(abund.prop)
  
  return(abund.prop)
} 




# unique plot function: get plots within both abund and ba data 1:N is the plot information
get.unique.plot <- function(data.abund, data.ba, N=8){
  data.abund <- as.data.frame(data.abund)
  data.ba    <- as.data.frame(data.ba)
  
  print('raw data');print(dim(data.abund));print(dim(data.ba))
  
  # select plots only with 1 row
  abund.plot <- names(table(data.abund$X))[which(table(data.abund$X)<2)]
  ba.plot    <- names(table(data.ba$X))[which(table(data.ba$X)<2)]
  print('only one row plots');print(length(abund.plot));print(length(ba.plot))
  
  # obtain unique plot id
  unique.plot <- abund.plot[which(abund.plot %in% ba.plot)]
  
  data.abund <- data.abund[which(data.abund$X %in% unique.plot),]
  data.ba    <- data.ba[which(data.ba$X %in% unique.plot),]
  print('final unique plots');print(dim(data.abund));print(dim(data.ba))
  
  # Standarize columns with 7
  sel.column <- c('X','Coords_y','Coords_x','Yr','Ps','B','D')
  data.abund <- cbind(data.abund[, sel.column], data.abund[, -c(1:N)])
  data.ba    <- cbind(data.ba[, sel.column],    data.ba[, -c(1:N)])
  
  output <- list(abund=data.abund, ba=data.ba)
  return(output)
}

# Fun 3: Match mycor & Get EcM proportion and diversity, 1:N is the plot information
get.pro.abund <- function(data1=fialand.abund, data2=fialand.ba, database='FIA', continent='North America'){
  
  data1 <- as.data.frame(data1)
  data2 <- as.data.frame(data2)
  
  ############# abundance data
  data.abund <- data1[, -c(1:7)]
  allsp <- data.frame(species.origin = colnames(data.abund), species= colnames(data.abund))
  
  # correct species name
  allsp$species <- gsub("\\."," ",  allsp$species)
  allsp$species <- gsub("_",  " ",  allsp$species)
  allsp$species <- gsub("-",  " ",  allsp$species)
  
  # TNRS name
  allsp <- merge(allsp, splist.TNRS[,1:8], by.x='species',by.y='species_submit',all.x=T, sort=F)
  allsp <- match.mycor(allsp)
  
  # test all species within splist.TNRS data
  print(paste0('all sp in data is ', nrow(allsp)))
  print(paste0('matched sp in TNRS is ', length(na.omit(allsp$species))))
  
  AM.sp <- subset(allsp, mycor=='AM')$species.origin
  EcM.sp <- subset(allsp, mycor=='EcM')$species.origin
  EcMAM.sp <- subset(allsp, mycor=='EcM_AM')$species.origin
  
  output <- data.frame(data1[,1:7],  
    AM.stem    = rowSums(data.abund[, AM.sp])+rowSums(data.abund[, EcMAM.sp])/2,
    EcM.stem   = rowSums(data.abund[, EcM.sp])+rowSums(data.abund[, EcMAM.sp])/2,
    EAM.stem   = rowSums(data.abund[, EcMAM.sp]),
    Total.stem = rowSums(data.abund),
    
    AM.SR   = specnumber(data.abund[, AM.sp])+specnumber(data.abund[, EcMAM.sp])/2,
    EcM.SR  = specnumber(data.abund[, EcM.sp])+specnumber(data.abund[, EcMAM.sp])/2,
    EAM.SR  = specnumber(data.abund[, EcMAM.sp]),
    SR      = specnumber(data.abund),
    shannon = diversity(data.abund, "shannon"),
    simpson = diversity(data.abund, "simpson"))
  
  output$AM.p.SR    <- with(output, AM.SR/SR)
  output$EcM.p.SR   <- with(output, EcM.SR/SR)
  output$EAM.p.SR   <- with(output, EAM.SR/SR)
  
  output$AM.p.stem  <- with(output, AM.stem/Total.stem)
  output$EcM.p.stem <- with(output, EcM.stem/Total.stem)
  output$EAM.p.stem <- with(output, EAM.stem/Total.stem)
  output$database   <- database
  output$continent  <- continent
  
  #############  BA data
  data.ba <- data2[, -c(1:7)]
  allsp <- data.frame(species.origin = colnames(data.ba), species= colnames(data.abund))
  
  # correct species name
  allsp$species <- gsub("\\."," ",  allsp$species)
  allsp$species <- gsub("_",  " ",  allsp$species)
  allsp$species <- gsub("-",  " ",  allsp$species)
  
  # TNRS name
  allsp <- merge(allsp, splist.TNRS[,1:8], by.x='species',by.y='species_submit',all.x=T, sort=F)
  allsp <- match.mycor(allsp)
  
  AM.sp <- subset(allsp, mycor=='AM')$species.origin
  EcM.sp <- subset(allsp, mycor=='EcM')$species.origin
  EcMAM.sp <- subset(allsp, mycor=='EcM_AM')$species.origin
  
  output.BA <- data.frame(X = data2[,'X'],  
    AM.ba    = rowSums(data.ba[, AM.sp])+rowSums(data.ba[, EcMAM.sp])/2,
    EcM.ba   = rowSums(data.ba[, EcM.sp])+rowSums(data.ba[, EcMAM.sp])/2,
    EAM.ba   = rowSums(data.ba[, EcMAM.sp]),
    Total.ba = rowSums(data.ba))
  
  output.BA$AM.p.BA    <- with(output.BA, AM.ba/Total.ba)
  output.BA$EcM.p.BA   <- with(output.BA, EcM.ba/Total.ba)
  output.BA$EAM.p.BA   <- with(output.BA, EAM.ba/Total.ba)
  output <- merge(output, output.BA, by='X', all.x=T, sort=F)
  
  output <- output[, c('X','Ps','Coords_y','Coords_x','D', 'SR','shannon','simpson',
    'AM.p.SR', 'AM.p.stem','AM.p.BA',
    'EcM.p.SR','EcM.p.stem','EcM.p.BA',
    'EAM.p.SR','EAM.p.stem','EAM.p.BA',
    'database','continent')]
  colnames(output)[1:5] <- c('Plot','Area','Lat','Lon','DBH')
  output <- na.omit(output)
  
  output$Plot <- paste0(database, '_', output$Plot)
  
  return(output)
  
}

# Fun 4: Read csv data and filter plots with recent year
read.ab <- function(file.name, encod="unknown"){
  path <- 'data/plot.data/'
  data <- data.table::fread(paste0(path, file.name), encoding = encod)
  names(data)[1] <- 'X'
  
  N1=nrow(data) # original rows
  
  data.plot <- data[,c('X','Yr')]
  
  pattern <- paste0('_', unique(data.plot$Yr))
  data.plot$plot <- stringi::stri_replace_all_regex(data.plot$X, pattern=pattern, replacement=c(''), vectorize=FALSE)
  data.plot <- data.plot %>% group_by(plot) %>% summarise(max.year = max(Yr,na.rm = T))
  plot_yr <- with(data.plot, paste0(plot, '_', max.year))
  
  output <- data[which(data$X %in% plot_yr), ]
  N2=nrow(output) # select plot number
  
  print(paste0(N1, ' plots to ', N2, ' plots'))
  
  return(output)
  
}

# Fun 5: Select plots randomly
sel.plot <- function(data, N){
  set.seed(123)
  data <- na.omit(data)
  output <- data[sample(1:nrow(data),N),]
  
  print(paste0(nrow(output), ' plots'))
  return(output)
}

# Fun 6: Down-sample data with 2*2 grid
sel.plot.grid <- function(data, N=3, cell=2){ # N is the max plots per grid; bin is the grid size
  data$latbin <- cut(data$Lat, seq(-90, 90, by=cell),  labels = seq(-90, 90, by=cell)[-1])
  data$lonbin <- cut(data$Lon, seq(-180, 180, by=cell),labels = seq(-180, 180, by=cell)[-1])
  data$grid   <- with(data, paste0(latbin, ".", lonbin))
  n.grid      <- length(unique(data$grid)) # number of grids
  print(paste0(n.grid, ' grids'))
  
  output <- data.frame()
  max.number <- N # set maximum plot number within a grid
  
  set.seed(123)
  for (i in 1:n.grid) {
    grid.i <- unique(data$grid)[i]
    data.i <- subset(data, grid==grid.i)
    
    if(nrow(data.i)>=max.number){
      out.i <- data.i[sample(1:nrow(data.i), size = max.number), ] # sample the maximum plot number
    }else{
      out.i <- data.i
    }
    
    output <- rbind(output, out.i)
  }
  
  print(paste0(nrow(output), ' plots'))
  return(output)
}

# abund data to long data
abund.long <- function(data, file){
  data <- as.data.frame(data)
  data <- melt(setDT(data), id.vars = 1:7, variable.name = "species")
  data <- as.data.frame(data)
  data <- subset(data, value > 0)
  data <- data[, c('X','species','Coords_y', 'Coords_x','Ps','value')]
  data <- merge(data, splist.TNRS[,1:8], by.x='species',by.y='species_submit',all.x=T, sort=F)
  data <- data[,c('X','speciesTNRS','genusTNRS','familyTNRS','Coords_y', 'Coords_x','Ps','value')]
  colnames(data) <- c('plot','speciesTNRS','genusTNRS','familyTNRS','Lat', 'Lon','area','density')
  data$plot <- paste0(file, data$plot)
  data$density_N <- with(data, area*density)
  #data <- unique(data)
  return(data)
}

# SEM function
sem.sr <- function(data=data, y.name='EcM.p.SR'){
  
  all.var <- c('SR', 
    'EcM.p.SR','EcM.p.stem','EcM.p.BA', 
    'bio1','AI')
  
  data <- data[,all.var] %>% na.omit() 
  
  data.use <- data[, all.var[-1]] %>% 
    #scale() %>% 
    as.data.frame() %>% # scale all predictors
    mutate(SR = data$SR)
  
  
  data.use$EcM  <- data.use[,y.name]
  data.use$EcM2 <- data.use$EcM^2
  
  mod <- psem(
    glm.nb(SR ~EcM+EcM2+bio1+AI, data = data.use),
    glm(EcM2~EcM       +bio1+AI, family=binomial(),data = data.use),
    glm(EcM~           +bio1+AI, family=binomial(),data = data.use), data=data.use  )
  
  
  mod.coef <- coefs(mod)
  mod.coef$Standcoef <- mod.coef$Std.Estimate
  
  mod.coef$Standcoef[1:4] <- with(data.use, c(
    mod.coef[1,3]* (sd(EcM)/sd(log(SR))),
    mod.coef[2,3]* (sd(EcM2)/sd(log(SR))),
    mod.coef[3,3]* (sd(bio1)/sd(log(SR))),
    mod.coef[4,3]* (sd(AI)/sd(log(SR)))
  ))
  
  mod.coef$Standcoef <- round(as.numeric(mod.coef$Standcoef), 4)
  print(mod.coef)
  
}
