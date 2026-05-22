######## R CODE: Behavioural thermoregulation and microclimate reshape climate-driven range forecasts #############
###################### Juan G. Rubalcaba*, Guillermo Fandos, José A. Díaz #########################################
###################### *jg.rubalcaba@gmail.com  ###################################################################

##### MICROCLIMATES and Te ----
# install.packages("C:/Users/juanv/Dropbox/Talento CAM/Field metabolic rates/rgdal_1.6-7.tar.gz", repos = NULL, type="source")
# require(devtools)
# install_github("ilyamaclean/microclima")

require(NicheMapR)
require(microclima)
require(stringr)
require(lubridate)

lat <- 40.523763
lon <- -3.768889
dstart <- "01/04/2022"
dfinish <- "31/07/2022"
Usrhyt <- 0.01

minshade <- 0
maxshade <- 75
Thcond <- 2.5
SpecHeat <- 870
Density <- 2.56
BulkDensity <- 1.3
windfac <- 1
REFL <- 0.2
cap <- FALSE
SLE <- 0.95
warm <- 0

clearsky <- FALSE

cat('downloading DEM via package elevatr /n')
dem <- microclima::get_dem(r = NA, lat = lat, lon = lon, resolution = 30, zmin = -20, xdims = 100, ydims = 100)
if(FALSE){
  elev <- raster::extract(dem, c(lon, lat))[1]
  xy <- data.frame(x = lon, y = lat)
  sp::coordinates(xy) = ~x + y
  sp::proj4string(xy) = "+init=epsg:4326"
  xy <- as.data.frame(sp::spTransform(xy, raster::crs(dem)))
  slope <- raster::terrain(dem, unit = "degrees")
  slope <- raster::extract(slope, xy)
  aspect <- raster::terrain(dem, opt = "aspect", unit = "degrees")
  aspect <- raster::extract(aspect, xy)
  ha36 <- 0
  for (i in 0:35) {
    har <- microclima::horizonangle(dem, i * 10, raster::res(dem)[1])
    ha36[i + 1] <- atan(raster::extract(har, xy)) * (180/pi)
  }
  hori <- spline(x = ha36, n = 24, method =  'periodic')$y
  hori[hori < 0] <- 0
  hori[hori > 90] <- 90
}else{
  slope <- 0
  aspect <- 0
  hori<- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
  ha36 <- spline(x = hori, n = 36, method =  'periodic')$y
  ha36[ha36 < 0] <- 0
  ha36[ha36 > 90] <- 90
  elev <- NA
}
soilgrids <- 0
spatial <- NA
ERR <- 1.5
micro <- micro_ncep(SLE = SLE, warm = warm, soilgrids = soilgrids, dstart = dstart, dfinish = dfinish,
                    Usrhyt = Usrhyt, slope = slope, aspect = aspect, REFL = REFL,
                    hori = hori, minshade = minshade, maxshade = maxshade,
                    loc = c(lon, lat), runshade = 1, run.gads = 1, snowmodel = 1,
                    BulkDensity =  BulkDensity, cap = cap,
                    Density = Density, Thcond = Thcond, SpecHeat = SpecHeat,
                    windfac = windfac, spatial = spatial, ERR = ERR, dem = dem)

soil <- as.data.frame(micro$soil)
metout <- as.data.frame(micro$metout)
shadsoil <- as.data.frame(micro$shadsoil)
shadmet <- as.data.frame(micro$shadmet)

require(lutz)
tz <- tz_lookup_coords(lat, lon)
dates <- as.POSIXct(as.character(micro$dates), format = "%Y-%m-%d %H:%M:%S", tz = tz)
dates2 <- as.POSIXct(as.character(micro$dates2), format = "%Y-%m-%d", tz = tz)

seas <- 'FALSE'
if(seas){
  dstart <- paste0("01/01/", '1987')
  dfinish <- paste0("30/06/", as.numeric('1987') + 1)
  dates2plot <- seq(as.POSIXct(paste0("01/07/", '1987',"00:00"), format = "%d/%m/%Y %H:%M", tz = tz), as.POSIXct(paste0("30/06/", as.numeric('1987') + 1,"23:00"), format = "%d/%m/%Y %H:%M", tz = tz), 3600)
  dates2plot2 <- seq(as.POSIXct(paste0("01/07/", '1987',"00:00"), format = "%d/%m/%Y %H:%M", tz = tz), as.POSIXct(paste0("30/06/", as.numeric('1987') + 1,"23:00"), format = "%d/%m/%Y %H:%M", tz = tz), 3600 * 24)
  if(as.POSIXct(dfinish, format = "%d/%m/%Y", tz = tz) > Sys.time()){
    dstart <- paste0("01/01/", '1987' - 1)
    dfinish <- paste0("30/06/", as.numeric('1987'))
    if(as.POSIXct(dfinish, format = "%d/%m/%Y", tz = tz) > Sys.time()){
      dfinish <- format(Sys.time() - 3600 * 48, "%d/%m/%Y", tz = tz)
    }
    dates2plot <- seq(as.POSIXct(paste0("01/07/", as.numeric('1987') - 1,"00:00"), format = "%d/%m/%Y %H:%M", tz = tz), as.POSIXct(paste0(dfinish,"23:00"), format = "%d/%m/%Y %H:%M", tz = tz), 3600)
    dates2plot2 <- seq(as.POSIXct(paste0("01/07/", as.numeric('1987') - 1,"00:00"), format = "%d/%m/%Y %H:%M", tz = tz), as.POSIXct(paste0(dfinish,"23:00"), format = "%d/%m/%Y %H:%M", tz = tz), 3600 * 24)
  }
}else{
  dstart <- paste0("01/07/", as.numeric('1987') - 1)
  dfinish <- paste0("31/12/", '1987')
  if(as.POSIXct(dfinish, format = "%d/%m/%Y") > Sys.time()){
    dfinish <- format(Sys.time() - 3600 * 48, "%d/%m/%Y", tz = tz)
  }
  dates2plot <- seq(as.POSIXct(paste0("01/01/", '1987', "00:00"), format = "%d/%m/%Y %H:%M", tz = tz), as.POSIXct(paste0(dfinish, "23:00"), format = "%d/%m/%Y %H:%M", tz = tz), 3600)
  dates2plot2 <- seq(as.POSIXct(paste0("01/01/", '1987',"00:00"), format = "%d/%m/%Y %H:%M", tz = tz), as.POSIXct(paste0(dfinish, "23:00"), format = "%d/%m/%Y %H:%M", tz = tz), 3600 * 24)
}
dates2plot2 <- as.POSIXct(format(dates2plot2, "%Y-%m-%d"), format = "%Y-%m-%d", tz = tz) # need this to avoid issue with daylight savings

metout$dates <- as.POSIXct(format(dates, "%Y/%m/%d"), format = "%Y/%m/%d", tz = tz)
soil$dates <- as.POSIXct(format(dates, "%Y/%m/%d"), format = "%Y/%m/%d", tz = tz)
shadmet$dates <- as.POSIXct(format(dates, "%Y/%m/%d"), format = "%Y/%m/%d", tz = tz)
shadsoil$dates <- as.POSIXct(format(dates, "%Y/%m/%d"), format = "%Y/%m/%d", tz = tz)

years <- as.numeric(format(dates, "%Y"))
daysinyear <- years
leapyears <- seq(1904, 3000, 4)
daysinyear[years %in% leapyears] <- 366
daysinyear[!years %in% leapyears] <- 365

write.csv(metout, file=paste0("C:/.../","ElPardo1997_metout.csv"))
write.csv(soil, file=paste0("C:/.../","ElPardo1997_soil.csv"))
write.csv(shadmet, file=paste0("C:/.../","ElPardo2022_shadmet.csv"))
write.csv(shadsoil, file=paste0("C:/.../","ElPardo2022_shadsoil.csv"))
