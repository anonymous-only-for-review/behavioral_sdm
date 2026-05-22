######## R CODE: Behavioural thermoregulation and microclimate reshape climate-driven range forecasts #############
###################### Juan G. Rubalcaba*, Guillermo Fandos, José A. Díaz #########################################
###################### *jg.rubalcaba@gmail.com  ###################################################################

####### Biophysical model functions ----
# Transient body temperature model
Tb_t_model <- function(M, # Body mass (g)
                       L, # Body length (m)
                       a, # Skin absorptivity (-)
                       ground=T, # Include conductive heat exchange to the ground (TRUE/FALSE)
                       use_Tsk=T,# Use sky temperature to model radiative heat exchange (TRUE/FALSE): if FALSE, uses air temperature
                       S,  # Solar radiation (W/m2)
                       Ta, # Air temperature (ºC)
                       Tsk,# Sky temperature (ºC)
                       Tg, # Air temperature (ºC)
                       v,  # Wind speed (m/s)
                       RH, # Relative humidity (%)
                       T0, # Initial body temperature (ºC)
                       time, # Time (s)
                       offset=0){ # Include temperature offset to simulate warming (increase in Ta)
  Ta <- Ta + offset
  
  # Derive morphological traits using allometric functions of body mass
  A = 0.0314 * pi * (M/1000)^(2/3)   # Surface area (m2) -> O'Connor 1999
  Ad = 0.4 * A # dorsal area (like in Fei et al. 2012 J Therm Biol)
  Ag = 0.6 * A # ventral area
  
  cp = 3.7 # Specific heat capacity body (J g-1 ºC-1) -> Porter et al. 1973
  C = M * cp
  
  ## Convective heat transfer coefficient
  k = 1.5207e-11*(Ta+273)^3 - 4.8574e-08*(Ta+273)^2 + 1.0184e-04*(Ta+273) - 3.9333e-04 # Thermal conductivity of air (W m-1 K-1)
  nu = -1.1555e-14*(Ta+273)^3 + 9.5728e-11*(Ta+273)^2 + 3.7604e-08*(Ta+273) - 3.4484e-06 # Kinematic viscosity of air (m2 s-1)
  Re = v * L / nu
  # Select one (Mitchell 1976)
  # Nu = 0.35 * Re^0.6 # Transverse to airflow
  # Nu = 0.1 * Re^0.74  # Parallel to airflow
  Nu = 1.36 * Re^0.39 # Postrate on surface, average parallel & perpendicular
  # Nu = 1.91 * Re^0.45 # Elevated from surface, average parallel & perpendicular
  hc = Nu * k / L # Convection heat transfer coef (W m-2 ?C-1)

  ## Conduction to the ground
  k_skin = 0.027 # Conductivity of the skin (W m-1 ?C-1)
  t_skin = 0.025 * (0.001 * M / (pi * 1000))^0.2 # thickness of skin (Stevenson 1985)
  hg = k_skin / t_skin # Conduction heat transfer coef (W m-2 ?C-1)
  
  if(!ground) hg = 0
  if(!use_Tsk) Tsk = Ta
  
  ## Thermal radiation
  epsilon = 0.98 # emissivity IR
  sigma = 5.67e-8 # Stefan Boltzmann constant (W m-2 K-4)
  Ra = 4 * epsilon * sigma * (Tsk+273)^3 # Radiative heat transfer coef (W m-2 ?C-1)
  Rg = 4 * epsilon * sigma * (Tg+273)^3
  
  ## Body temperature (see Supplementary Information 1; eqs.7-9)
  j = Ad / C * (a * S + Ra * Tsk + hc * Ta) + Ag / C * Tg * (Rg + hg)
  theta = Ad / C * (Ra + hc) + Ag / C * (Rg + hg)
  
  Te <- j/theta
  Tb <- Te + (T0 - Te) * exp(-theta * time)
  return(Tb)
}
# Operative temperature model
Te_model <- function(M, L, a, ground=T, use_Tsk=T, S, Ta, Tsk, Tg, v, RH, offset=0){
  Ta <- Ta + offset
  
  # Derive morphological traits using allometric functions of body mass
  A = 0.0314 * pi * (M/1000)^(2/3)   # Surface area (m2) -> O'Connor 1999
  Ad = 0.4 * A # dorsal area (like in Fei et al. 2012 J Therm Biol)
  Ag = 0.6 * A # ventral area
  
  cp = 3.7 # Specific heat capacity body (J g-1 ºC-1) -> Porter et al. 1973
  C = M * cp
  
  ## Convective heat transfer coefficient
  k = 1.5207e-11*(Ta+273)^3 - 4.8574e-08*(Ta+273)^2 + 1.0184e-04*(Ta+273) - 3.9333e-04 # Thermal conductivity of air (W m-1 K-1)
  nu = -1.1555e-14*(Ta+273)^3 + 9.5728e-11*(Ta+273)^2 + 3.7604e-08*(Ta+273) - 3.4484e-06 # Kinematic viscosity of air (m2 s-1)
  Re = v * L / nu
  # Select one (Mitchell 1976)
  # Nu = 0.35 * Re^0.6 # Transverse to airflow
  # Nu = 0.1 * Re^0.74  # Parallel to airflow
  Nu = 1.36 * Re^0.39 # Postrate on surface, average parallel & perpendicular
  # Nu = 1.91 * Re^0.45 # Elevated from surface, average parallel & perpendicular
  hc = Nu * k / L # Convection heat transfer coef (W m-2 ?C-1)
  
  ## Conduction to the ground
  k_skin = 0.027 # Conductivity of the skin (W m-1 ?C-1)
  t_skin = 0.025 * (0.001 * M / (pi * 1000))^0.2 # thickness of skin (Stevenson 1985)
  hg = k_skin / t_skin # Conduction heat transfer coef (W m-2 ?C-1)
  
  if(!ground) hg = 0
  if(!use_Tsk) Tsk = Ta
  
  ## Thermal radiation
  epsilon = 0.98 # emissivity IR
  sigma = 5.67e-8 # Stefan Boltzmann constant (W m-2 K-4)
  Ra = 4 * epsilon * sigma * (Tsk+273)^3 # Radiative heat transfer coef (W m-2 ?C-1)
  Rg = 4 * epsilon * sigma * (Tg+273)^3
  
  ## Body temperature (see Supplementary Information 1; eqs.7-9)
  j = Ad / C * (a * S + Ra * Tsk + hc * Ta) + Ag / C * Tg * (Rg + hg)
  theta = Ad / C * (Ra + hc) + Ag / C * (Rg + hg)
  
  return(j/theta)
}

####### Load microclimate and field temperature data ----

dir <- "C:/..."
year="2022" # Select year ("1997" OR "2022")
datamicro_soil <- read.csv(paste0(dir,"/Microclimate data/ElPardo",year,"_soil.csv"), sep=",", header=T) # FROM http://bioforecasts.science.unimelb.edu.au/app_direct/soil_ncep/
datamicro_metout <- read.csv(paste0(dir,"/Microclimate data/ElPardo",year,"_metout.csv"), sep=",", header=T)
datamicro_shadmet <- read.csv(paste0(dir,"/Microclimate data/ElPardo",year,"_shadmet.csv"), sep=",", header=T)
datamicro_shadsoil <- read.csv(paste0(dir,"/Microclimate data/ElPardo",year,"_shadsoil.csv"), sep=",", header=T)

datamicro_metout$day <- as.numeric(format(as.Date(datamicro_metout$dates,format="%Y-%m-%d"), format = "%d"))
datamicro_metout$month <- format(as.Date(datamicro_metout$dates,format="%Y-%m-%d"), format = "%m")
datamicro_metout$year <- format(as.Date(datamicro_metout$dates,format="%Y-%m-%d"), format = "%Y")
datamicro_metout$datetime <- as.POSIXct(paste(datamicro_metout$dates, datamicro_metout$TIME/1440*24), format="%Y-%m-%d %H")

datamicro_shadmet$day <- as.numeric(format(as.Date(datamicro_shadmet$dates,format="%Y-%m-%d"), format = "%d"))
datamicro_shadmet$month <- format(as.Date(datamicro_shadmet$dates,format="%Y-%m-%d"), format = "%m")
datamicro_shadmet$year <- format(as.Date(datamicro_shadmet$dates,format="%Y-%m-%d"), format = "%Y")
datamicro_shadmet$datetime <- as.POSIXct(paste(datamicro_shadmet$dates, datamicro_shadmet$TIME/1440*24), format="%Y-%m-%d %H")

datamicro_soil$day <- as.numeric(format(as.Date(datamicro_soil$dates,format="%Y-%m-%d"), format = "%d"))
datamicro_soil$month <- format(as.Date(datamicro_soil$dates,format="%Y-%m-%d"), format = "%m")
datamicro_soil$year <- format(as.Date(datamicro_soil$dates,format="%Y-%m-%d"), format = "%Y")
datamicro_soil$datetime <- as.POSIXct(paste(datamicro_soil$dates, datamicro_soil$TIME/1440*24), format="%Y-%m-%d %H")

datamicro_shadsoil$day <- as.numeric(format(as.Date(datamicro_shadsoil$dates,format="%Y-%m-%d"), format = "%d"))
datamicro_shadsoil$month <- format(as.Date(datamicro_shadsoil$dates,format="%Y-%m-%d"), format = "%m")
datamicro_shadsoil$year <- format(as.Date(datamicro_shadsoil$dates,format="%Y-%m-%d"), format = "%Y")
datamicro_shadsoil$datetime <- as.POSIXct(paste(datamicro_shadsoil$dates, datamicro_shadsoil$TIME/1440*24), format="%Y-%m-%d %H")

data_fieldTb <- read.csv(paste0(dir,"/Microclimate data/field_Tb.csv"), sep=",", header=T)
data_fieldTe <- read.csv(paste0(dir,"/Microclimate data/field_Te.csv"), sep=",", header=T)

####### Operative temperature model ----
# Lizard's traits
M = 8.3 # Body mass (g)
L = 0.069 # SVL (m)
a = 0.9  # Skin absorbance (P. muralis; Clusella-Trullas et al. 2011)

maxshade = 75
shlev = 1-maxshade/100 # % sun in shaded spots
datamicro_metout$Te_sun <- Te_model(M, L, a, use_Tsk = T, ground = F, S=datamicro_metout$SOLR, Ta=datamicro_metout$TALOC, 
                                    Tsk=datamicro_metout$TSKYC, Tg=datamicro_soil$D0cm, v=datamicro_metout$VLOC,
                                    RH=datamicro_metout$RHLOC)
datamicro_shadmet$Te_shade <- Te_model(M, L, a, use_Tsk = T, ground = F, S=datamicro_shadmet$SOLR*shlev, Ta=datamicro_shadmet$TALOC, 
                                       Tsk=datamicro_shadmet$TSKYC, Tg=datamicro_shadsoil$D0cm, v=datamicro_shadmet$VLOC,
                                       RH=datamicro_shadmet$RHLOC)

####### Body temperature model ----

# Behavioral thermoregulation parameters
lambda = 0.8  # Thermoregulatory constraint (-)
Tpref = 33.8    # Preferred temperature (ºC) -> Tpref_MAY -> 32.9; 
lag = 60       # time lag (sec)
maxshade = 75 # maximum shade level
shlev = 1-maxshade/100 # % sun in shaded spots
# Year	Lambda	Tpref
# 1997	 0.95	  33.8
# 2022   0.8	  33.8

# Lizard's traits
M = 8.3 # Body mass (g)
L = 0.069 # SVL (m)
a = 0.9  # Skin absorbance (P. muralis; Clusella-Trullas et al. 2011)

days <- sort(unique(data_fieldTb$day[which(data_fieldTb$month==5 & data_fieldTb$year==year)])) # Days with field Tb information
reps <- 10 # number of simulations
Tbs <- locations <- array(NA, dim=c(1,reps))
for(day in days){
  # Subset of the microclimate database to consider only days in which empirical Tb data were collected
  datasun <- datamicro_metout[which(datamicro_metout$month=="05" & datamicro_metout$day == day),]
  datashade <- datamicro_shadmet[which(datamicro_shadmet$month=="05" & datamicro_shadmet$day == day),]
  datasoil <- datamicro_soil[which(datamicro_soil$month=="05" & datamicro_soil$day == day),]
  datasoilshade <- datamicro_shadsoil[which(datamicro_shadsoil$month=="05" & datamicro_shadsoil$day == day),]
  
  seq_min <- seq(1,nrow(datasun),length.out=nrow(datasun)*60)
  TIME_fun <- approxfun(datasun$TIME ~ I(1:nrow(datasun)))
  TIME <- TIME_fun(seq_min)
  
  # transform hourly estimations into minutely values to run the behavioral
  # thermoregulation model
  SOLR_sun_fun <- approxfun(datasun$SOLR ~ I(1:nrow(datasun)))
  SOLR_sun <- SOLR_sun_fun(seq_min)
  TA_sun_fun <- approxfun(datasun$TALOC ~ I(1:nrow(datasun)))
  TA_sun <- TA_sun_fun(seq_min)
  TS_sun_fun <- approxfun(datasoil$D0cm ~ I(1:nrow(datasoil)))
  TS_sun <- TS_sun_fun(seq_min)
  TSK_sun_fun <- approxfun(datasun$TSKYC ~ I(1:nrow(datasun)))
  TSK_sun <- TSK_sun_fun(seq_min)
  V_sun_fun <- approxfun(datasun$VLOC ~ I(1:nrow(datasun)))
  V_sun <- V_sun_fun(seq_min)
  
  TA_shade_fun <- approxfun(datashade$TALOC ~ I(1:nrow(datashade)))
  TA_shade <- TA_shade_fun(seq_min)
  TS_shade_fun <- approxfun(datasoilshade$D0cm ~ I(1:nrow(datasoilshade)))
  TS_shade <- TS_shade_fun(seq_min)
  TSK_shade_fun <- approxfun(datashade$TSKYC ~ I(1:nrow(datashade)))
  TSK_shade <- TSK_shade_fun(seq_min)
  V_shade_fun <- approxfun(datashade$VLOC ~ I(1:nrow(datashade)))
  V_shade <- V_shade_fun(seq_min)
  
  # Behavioral thermoregulation model
  duration <- 23*60 # min
  Tb <- location <- array(NA,dim=c(duration,reps))
  Tb[1,] <- TA_shade[1]
  location[1,] <- -1
  for(rep in 1:reps){ # Repeat "rep" times
    for(i in 2:duration){ # for each minute:
      # Derive expected Tb in both sun and shaded conditions
      TbSun <- Tb_t_model(M=M, L=L, a=a, ground=F, use_Tsk=T, S=SOLR_sun[i-1], Ta=TA_sun[i-1],
                          Tsk=TSK_sun[i-1],Tg=TS_sun[i-1], v=V_sun[i-1], RH=0.5,
                          T0=Tb[i-1,rep], time=lag)
      TbShade <- Tb_t_model(M=M, L=L, a=a, ground=F, use_Tsk=T, S=SOLR_sun[i-1]*shlev, Ta=TA_shade[i-1],
                            Tsk=TSK_shade[i-1],Tg=TS_shade[i-1], v=V_shade[i-1], RH=0.5,
                            T0=Tb[i-1,rep], time=lag)
      
      # Compute abs distance to Tpref in both sun and shaded conditions
      wTsun <- abs(TbSun - Tpref) 
      wTshade <- abs(TbShade - Tpref) 
      Z <- exp(-lambda*wTsun) + exp(-lambda*wTshade)
      if(Z==0) Z=1e-3
      Psun <- exp(-lambda*wTsun) / Z # probability of selecting sun-exposed conditions
      
      if(rbinom(1,1,Psun)){ # if the animal moves to sun
        Tb[i,rep] <- TbSun
        location[i,rep] <- 1
      }else{ # if it moves to the shade
        Tb[i,rep] <- TbShade
        location[i,rep] <- -1
      }
    }
  }
  Tbs <- rbind(Tbs, Tb)
  locations <- rbind(locations, location)
}
Tbs <- Tbs[-1,]
locations <- locations[-1,]
ndays <- sort(rep(days,duration))
ntimes <- rep(sort(rep(2:24, 60)),length(days))

days_field_data <- as.numeric(days)+120 # Calculate DOY
data_sub_sun <- datamicro_metout[which(datamicro_metout$DOY == days_field_data[1]),]
data_sub_shade <- datamicro_shadmet[which(datamicro_shadmet$DOY == days_field_data[1]),]
data_fieldTe_SUN <- data_fieldTe[which(data_fieldTe$month==5 & data_fieldTe$year == year & data_fieldTe$SUN_EXPOSURE == "Sun"),]
data_fieldTe_SUN$DOY <- as.numeric(data_fieldTe_SUN$day)+120
data_fieldTe_SUN_DOY <- data_fieldTe_SUN[which(data_fieldTe_SUN$DOY == days_field_data[1]),]
data_fieldTe_SHADE <- data_fieldTe[which(data_fieldTe$month==5 & data_fieldTe$year == year & data_fieldTe$SUN_EXPOSURE == "Shade"),]
data_fieldTe_SHADE$DOY <- as.numeric(data_fieldTe_SHADE$day)+120
data_fieldTe_SHADE_DOY <- data_fieldTe_SHADE[which(data_fieldTe_SHADE$DOY == days_field_data[1]),]
data_fieldTe_SUNSHADE <- data_fieldTe[which(data_fieldTe$month==5 & data_fieldTe$year == year & data_fieldTe$SUN_EXPOSURE == "Sun-shade"),]
data_fieldTe_SUNSHADE$DOY <- as.numeric(data_fieldTe_SUNSHADE$day)+120
data_fieldTe_SUNSHADE_DOY <- data_fieldTe_SUNSHADE[which(data_fieldTe_SUNSHADE$DOY == days_field_data[1]),]

for(i in 2:length(days_field_data)){
  data_sub_sun <- rbind(data_sub_sun,datamicro_metout[which(datamicro_metout$DOY == days_field_data[i]),])
  data_sub_shade <- rbind(data_sub_shade,datamicro_shadmet[which(datamicro_shadmet$DOY == days_field_data[i]),])
  data_fieldTe_SUN_DOY <- rbind(data_fieldTe_SUN_DOY,data_fieldTe_SUN[which(data_fieldTe_SUN$DOY == days_field_data[i]),])
  data_fieldTe_SHADE_DOY <- rbind(data_fieldTe_SHADE_DOY,data_fieldTe_SHADE[which(data_fieldTe_SHADE$DOY == days_field_data[i]),])
  data_fieldTe_SUNSHADE_DOY <- rbind(data_fieldTe_SUNSHADE_DOY,data_fieldTe_SUNSHADE[which(data_fieldTe_SUNSHADE$DOY == days_field_data[i]),])
}

shift = 2 # adjust summer time difference
code_modTe <- paste0(data_sub_sun$DOY,"_", shift+data_sub_sun$TIME/60)
code_obsTe_sun <- paste0(data_fieldTe_SUN_DOY$DOY,"_", data_fieldTe_SUN_DOY$TIME)
code_obsTe_shade <- paste0(data_fieldTe_SHADE_DOY$DOY,"_", data_fieldTe_SHADE_DOY$TIME)
code_obsTe_sunshade <- paste0(data_fieldTe_SUNSHADE_DOY$DOY,"_", data_fieldTe_SUNSHADE_DOY$TIME)

modTe_sun <- data_sub_sun$Te[match(code_obsTe_sun,code_modTe)]
modTe_shade <- data_sub_shade$Te[match(code_obsTe_shade,code_modTe)]
modTe_sunshade <- data_sub_shade$Te[match(code_obsTe_sunshade,code_modTe)]

# Summarize Observed vs predicted operative temperatures for each hour
Te_sun_obs <- tapply(data_fieldTe_SUN_DOY$Te, data_fieldTe_SUN_DOY$TIME, mean)
Te_sun_obs_SD <- tapply(data_fieldTe_SUN_DOY$Te, data_fieldTe_SUN_DOY$TIME, sd)
Te_shade_obs <- tapply(data_fieldTe_SHADE_DOY$Te, data_fieldTe_SHADE_DOY$TIME, mean)
Te_shade_obs_SD <- tapply(data_fieldTe_SHADE_DOY$Te, data_fieldTe_SHADE_DOY$TIME, sd)
Te_sun_pred <- tapply(modTe_sun, data_fieldTe_SUN_DOY$TIME, mean)
Te_sun_pred_SD <- tapply(modTe_sun, data_fieldTe_SUN_DOY$TIME, sd)
Te_shade_pred <- tapply(modTe_shade, data_fieldTe_SHADE_DOY$TIME, mean)
Te_shade_pred_SD <- tapply(modTe_shade, data_fieldTe_SHADE_DOY$TIME, sd)
n_obs_sun <- tapply(data_fieldTe_SUN_DOY$Te, data_fieldTe_SUN_DOY$TIME, length)
n_obs_shade <- tapply(data_fieldTe_SHADE_DOY$Te, data_fieldTe_SHADE_DOY$TIME, length)

# Operative temperatures in the sun (Observed vs Predicted)
ggplot(mapping=aes(y=Te_sun_obs, x=Te_sun_pred)) + # 500 x 400
  theme_classic() + ylab("Observed Te (ºC)") + xlab("Predicted Te (ºC)") +
  ylim(20,80) + xlim(20,80) +
  theme(axis.text = element_text(size=10, colour="black"),
        axis.title = element_text(size=12)) +
  geom_point(size=n_obs_sun*0.1) + geom_abline(intercept=0, slope = 1) +
  geom_errorbar(aes(ymin=Te_sun_obs-Te_sun_obs_SD, ymax=Te_sun_obs+Te_sun_obs_SD), width=1, position=position_dodge(.9)) +
  geom_errorbar(aes(xmin=Te_sun_pred-Te_sun_pred_SD, xmax=Te_sun_pred+Te_sun_pred_SD), width=1, position=position_dodge(.9)) 

# Operative temperatures in the shade (Observed vs Predicted)
ggplot(mapping=aes(y=Te_shade_obs, x=Te_shade_pred)) + # 500 x 400
  theme_classic() + ylab("Observed Te (ºC)") + xlab("Predicted Te (ºC)") +
  ylim(10,40) + xlim(10,40) +
  theme(axis.text = element_text(size=10, colour="black"),
        axis.title = element_text(size=12)) +
  geom_point(size=n_obs_shade*0.1) + geom_abline(intercept=0, slope = 1) +
  geom_errorbar(aes(ymin=Te_shade_obs-Te_shade_obs_SD, ymax=Te_shade_obs+Te_shade_obs_SD), width=1, position=position_dodge(.9)) +
  geom_errorbar(aes(xmin=Te_shade_pred-Te_shade_pred_SD, xmax=Te_shade_pred+Te_shade_pred_SD), width=1, position=position_dodge(.9)) 

# Summarize Observed vs predicted body temperatures for each hour
predTb <- tapply(Tbs[,1], ntimes, function(x) mean(x))
predTb_SD <- tapply(Tbs[,1], ntimes, function(x) sd(x))
obsTb <- tapply(data_fieldTb$Tb, data_fieldTb$TIME, function(x) mean(x,na.rm=T))
obsTb_SD <- tapply(data_fieldTb$Tb, data_fieldTb$TIME, function(x) sd(x,na.rm=T))
n_obs <- tapply(data_fieldTb$Tb, data_fieldTb$TIME, length)

y <- predTb[match(names(obsTb), names(predTb))]
predTb <- y[which(!is.na(y))]
y <- predTb_SD[match(names(obsTb), names(predTb_SD))]
predTb_SD <- y[which(!is.na(y))]

ggplot(mapping=aes(y=obsTb, x=predTb)) +
  theme_classic() + ylab("Observed Tb (ºC)") + xlab("Predicted Tb (ºC)") +
  ylim(20,40) + xlim(20,40) +
  theme(axis.text = element_text(size=10, colour="black"),
        axis.title = element_text(size=12)) +
  geom_point(size=n_obs*0.5) + geom_abline(intercept=0, slope = 1) +
  geom_errorbar(aes(ymin=obsTb-obsTb_SD, ymax=obsTb+obsTb_SD), width=.2, position=position_dodge(.9)) +
  geom_errorbar(aes(xmin=predTb-predTb_SD, xmax=predTb+predTb_SD), width=.2, position=position_dodge(.9)) 


####### Geographical projections ----
require(NicheMapR)
require(microclima)
require(raster)
require(tidyr)
require(dplyr)

dir <- "C:/..." # Set directory
load(paste0(dir,"Sources/xy.values.RData"))
load(paste0(dir,"Sources/cells.RData"))
map <- raster(paste0(dir,"Sources/map.grd"))

# Create a function to run the NicheMapR microcliamte model across cells
# extract conditions in the sun and in the shade
# transform hourly estimations into minutely values to run the behavioral
# thermoregulation model
compute_microclimates <- function(lat, long, month, maxshade, warm=0){
  ## Run NicheMapR microcolimate model (Kearney et al. 2017 Ecography)
  xy.values <- c(x=long, y=lat)
  
  micro <- micro_global(run.gads=0, loc = xy.values, minshade = minshade, maxshade = maxshade, Usrhyt = 0.01, warm = warm)
  
  microclim_sun <- as.data.frame(micro$metout)    # microclimatic conditions in the sun
  microclim_shade <- as.data.frame(micro$shadmet) # and in the shade
  
  soil_sun <- as.data.frame(micro$soil)          # soil temperature in the sun
  soil_shade <- as.data.frame(micro$shadsoil)    # and in the shade
  
  ## Extract variables
  
  Ta_sun = microclim_sun$TALOC # Air temperature (?C)
  Ta_shade = microclim_shade$TALOC
  
  Tg_sun = soil_sun$D0cm      # Soil surface temperature
  Tg_shade = soil_shade$D0cm
  
  S_sun = microclim_sun$SOLR  # Solar radiation (Wm-2)
  S_shade = microclim_shade$SOLR * (1 - maxshade/100) 
  
  V_sun = microclim_sun$VLOC # Wind velocity (ms-2)
  V_shade = microclim_shade$VLOC
  
  Tsky_sun = microclim_sun$TSKYC # Sky temperature (thermal radiation)
  Tsky_shade = microclim_shade$TSKYC
  
  DOY = microclim_sun$DOY # day of the year (we will use it later)
  ZEN = microclim_sun$ZEN # zenith angle
  
  microclimate_sun_df <- data.frame(DOY=DOY, ZEN=ZEN, Ta=Ta_sun, Tg=Tg_sun, S=S_sun, V=V_sun, Tsky=Tsky_sun)
  microclimate_shade_df <- data.frame(DOY=DOY, ZEN=ZEN, Ta=Ta_shade, Tg=Tg_shade, S=S_shade, V=V_shade, Tsky=Tsky_shade)
  
  # filter dataframes by the month of interest
  
  Jday = c(15, 46, 74, 105, 135, 166, 196, 227, 258, 288, 319, 349)[month] # transform into Julian day
  
  microclimate_sun_df = microclimate_sun_df %>% filter(DOY == Jday) # filter by Julian day
  microclimate_shade_df = microclimate_shade_df %>% filter(DOY == Jday)
  
  # Temporal interpolation (increase resolution: hours -> minutes)
  
  rows <- nrow(microclimate_sun_df)
  
  microclimate_sun_min <- data.frame(array(NA, dim=c(60 * rows, 7)))
  for(i in 2:7){
    var_hour = microclimate_sun_df[,i]
    var_hour_interp_fun = approxfun(var_hour ~ 1:rows)
    microclimate_sun_min[,i] <- var_hour_interp_fun(seq(1, rows, length.out = 60*rows))
  }
  microclimate_shade_min <- data.frame(array(NA, dim=c(60 * rows, 7)))
  for(i in 2:7){
    var_hour = microclimate_shade_df[,i]
    var_hour_interp_fun = approxfun(var_hour ~ 1:rows)
    microclimate_shade_min[,i] <- var_hour_interp_fun(seq(1, rows, length.out = 60*rows))
  }
  microclimate_sun_min[,1] <- microclimate_shade_min[,1] <- sort(rep(microclimate_sun_df$DOY, 60))
  colnames(microclimate_sun_min) <- colnames(microclimate_shade_min) <- colnames(microclimate_sun_df)
  
  return(list(microclimate_sun_min=microclimate_sun_min, microclimate_shade_min=microclimate_shade_min))
}

# Load climate change projection

load(file=paste0(dir,"Sources/predicted_warming.RData")) # MIROC6, RCP 8.5, 2041-2060
str(predicted_warming) # col1: predicted temperature increase in may
                       # col2: predicted T increase in june

####### Behavioral thermoregulation model ----

# General parameters
M = 8.3 # Body mass (g)
L = 0.069 # SVL (m)
a = 0.9  # Skin absorbance (P. muralis; Clusella-Trullas et al. 2011)
lag = 60       
reps = 10     
maxshade = 75
minshade = 0

# These parameters need to be adjusted to either may or june
month = 6
lambda = 1.4 # Thermoregulation parameter (may: 0.9; june: 1.4)   
Tpref = 35.7 # Preferred temperature (may: 33.8, june: 35.7)
warm = predicted_warming[,2] 

meanTb <- daytimeTb <- meanTe_sun <- meanTe_shade <- daytimeTe_sun <- daytimeTe_shade <- E_daytime <- db_daytime <- Activity <- numeric(nrow(xy.values))
for(cell in 1:nrow(xy.values)){
  warm_cell <- warm[cell] # To simulate current conditions (no warming), set: warm_cell <- 0
  
  micro <- compute_microclimates(lat=xy.values[cell,2], long=xy.values[cell,1], month=month, maxshade=maxshade, warm=warm_cell)
  datasun <- micro$microclimate_sun_min
  datashade <- micro$microclimate_shade_min
  duration = nrow(datasun)
  
  TeSun <- Te_model(M=M, L=L, a=a, ground=F, use_Tsk=T, S=datasun$S, Ta=datasun$Ta,
                    Tsk=datasun$Tsky,Tg=datasun$Tg, v=datasun$V, RH=0.5)
  TeShade <- Te_model(M=M, L=L, a=a, ground=F, use_Tsk=T, S=datashade$S, Ta=datashade$Ta,
                      Tsk=datashade$Tsky,Tg=datashade$Tg, v=datashade$V, RH=0.5)
  max(TeSun)
  
  Tb <- location <- array(NA,dim=c(24*60,reps))
  Tb[1,] <- TeShade[1]
  location[1,] <- -1
  for(rep in 1:reps){
    for(i in 2:duration){
      TbSun <- Tb_t_model(M=M, L=L, a=a, ground=F, use_Tsk=T, S=datasun$S[i-1], Ta=datasun$Ta[i-1],
                          Tsk=datasun$Tsky[i-1],Tg=datasun$Tg[i-1], v=datasun$V[i-1], RH=0.5,
                          T0=Tb[i-1,rep], time=lag)
      TbShade <- Tb_t_model(M=M, L=L, a=a, ground=F, use_Tsk=T, S=datashade$S[i-1], Ta=datashade$Ta[i-1],
                            Tsk=datashade$Tsky[i-1],Tg=datashade$Tg[i-1], v=datashade$V[i-1], RH=0.5,
                            T0=Tb[i-1,rep], time=lag)
      
      # curve(1-dnorm(x, mean=32.9, sd=1.476),10,40)
      
      wTsun <- abs(TbSun - Tpref) # 1-dnorm(TbSun, mean=Tpref, sd=1) 
      wTshade <- abs(TbShade - Tpref) #1-dnorm(TbShade, mean=32.9, sd=1.476)
      Z <- exp(-lambda*wTsun) + exp(-lambda*wTshade)
      if(Z==0) Z=1e-3
      Psun <- exp(-lambda*wTsun) / Z
      
      if(rbinom(1,1,Psun)){
        Tb[i,rep] <- TbSun
        location[i,rep] <- 1
      }else{
        Tb[i,rep] <- TbShade
        location[i,rep] <- -1
      }
    }
  }
  
  # Mean daytime Tb and Te
  ## thermoregulation window
  x <- which(TeSun>Tpref)
  if(length(x)==0){
    sunrise <- which.max(TeSun)
    sunset <- which.max(TeSun)
  } else {
    sunrise <- min(x) 
    sunset <- max(x) 
  }
  
  meanTb[cell] <- mean(Tb)
  meanTb_time <- tapply(rowMeans(Tb), 1:nrow(datasun), mean)
  daytimeTb[cell] <- mean(meanTb_time[sunrise:sunset])
  
  meanTe_sun[cell] <- mean(TeSun)
  meanTe_shade[cell] <- mean(TeShade)
  meanTe_sun_time <- tapply(TeSun, 1:nrow(datasun), mean)
  meanTe_shade_time <- tapply(TeShade, 1:nrow(datasun), mean)
  daytimeTe_sun[cell] <- mean(meanTe_sun_time[sunrise:sunset])
  daytimeTe_shade[cell] <- mean(meanTe_shade_time[sunrise:sunset])
  
  # Thermoregulatory effectiveness (Hertz et al. 1993)
  db <- rowMeans(abs(Tb - Tpref)) # thermoregulatory inaccuracy
  meanTe <- rowMeans(cbind(TeSun,TeShade))
  de <- rowMeans(abs(cbind(TeSun,TeShade)-Tpref))
  
  E_total <- 1 - db/de
  
  db_daytime[cell] <- mean(db[sunrise:sunset])
  E_daytime[cell] <- mean(E_total[sunrise:sunset])
  
  # Activity (minutes for active thermoregulation)
  Activity[cell] <- sunset-sunrise
  
  print(cell/nrow(xy.values)*100)
}

# Store may data
dataMAY_warm <- data.frame(meanTb, daytimeTb, meanTe_sun, meanTe_shade, daytimeTe_sun,
                       daytimeTe_shade, E_daytime, db_daytime, Activity, xy.values)
# Store june data
dataJUNE_warm <- data.frame(meanTb, daytimeTb, meanTe_sun, meanTe_shade, daytimeTe_sun,
                           daytimeTe_shade, E_daytime, db_daytime, Activity, xy.values)

####### Load model output ----
# Simulated temperatures for May under current conditions
dataMAY <- read.table(file=paste0(dir,"Sources/dataMAY.R"))
# May temperatures under climate change scenario 
dataMAY_warm <- read.table(file=paste0(dir,"Sources/dataMAY_warm.R"))
# June temperatures current conditions
dataJUNE <- read.table(file=paste0(dir,"Sources/dataJUNE.R"))
# June temperatures climate change scenario 
dataJUNE_warm <- read.table(file=paste0(dir,"Sources/dataJUNE_warm.R"))

## Simulation data processing
# Mean thermoregulatory performance (temperature innacuracy )
meandb <- rowMeans(cbind(dataMAY$db_daytime, dataJUNE$db_daytime)) # current 
meandb_warm <- rowMeans(cbind(dataMAY_warm$db_daytime, dataJUNE_warm$db_daytime)) # future

# Mean thermoregulatory window
meanActivity <- rowMeans(cbind(dataMAY$Activity, dataJUNE$Activity)) # current 
meanActivity_warm <- rowMeans(cbind(dataMAY_warm$Activity, dataJUNE_warm$Activity)) # future 

# Predicted temperature increase in May and June
TbIncrease_MAY <- dataMAY_warm$meanTb - dataMAY$meanTb # Body temperature
TeIncrease_MAY <- dataMAY_warm$meanTe_shade - dataMAY$meanTe_shade # Operative temperature
TbIncrease_JUNE <- dataJUNE_warm$meanTb - dataJUNE$meanTb 
TeIncrease_JUNE <- dataJUNE_warm$meanTe_shade - dataJUNE$meanTe_shade 
bufferMAY <- TbIncrease_MAY / TeIncrease_MAY # Predicted buffering
bufferJUNE <- TbIncrease_JUNE / TeIncrease_JUNE

####### Generate raster layers ----

buffer_map <- meandb_map <- meandb_warm_map <- meanActivity_map <- meanActivity_warm_map <- map
for(i in 1:length(cells)){
  buffer_map[cells[i]] <- mean(c(bufferMAY[i], bufferJUNE[i]))
  meandb_map[cells[i]] <- meandb[i]
  meandb_warm_map[cells[i]] <- meandb_warm[i] 
  meanActivity_map[cells[i]] <- meanActivity[i]
  meanActivity_warm_map[cells[i]] <- meanActivity_warm[i]
}

writeRaster(buffer_map, file=paste0(dir,"/Sources/buffer_map.grd"), overwrite=T)
writeRaster(meandb_map, file=paste0(dir,"/Sources/meandb_map.grd"), overwrite=T)
writeRaster(meandb_warm_map, file=paste0(dir,"/Sources/meandb_warm_map.grd"), overwrite=T)
writeRaster(meanActivity_map, file=paste0(dir,"/Sources/meanActivity_map.grd"), overwrite=T)
writeRaster(meanActivity_warm_map, file=paste0(dir,"/Sources/meanActivity_warm_map.grd"), overwrite=T)
