### MB2 final project Daniel Gruschwitz
### (semi) automtatic species distribution modelling
### Environment setup #########################################################

rp <- c("raster", "sf", "RColorBrewer", "rgeos", "Rsagacmd","randomForest", "verification", "dplyr", "rfUtilities", "corrplot", "stringr", "wellknown", "tidyverse","rgbif", "devtools",  "sp", "rgdal") # required packages
ip <- rp[!(rp %in% installed.packages()[, "Package"])] # subset packages

if(length(ip)) install.packages(ip, dependencies = TRUE) # install packages
lapply(rp, require, character.only = TRUE) # load packages

rm(rp, ip) # remove vector

rp <- c("gdalUtils", "scrubr")
rp2 <- c("gearslaboratory/gdalUtils", "ropensci/scrubr")
ip <- rp[!(rp %in% installed.packages()[, "Package"])] # subset packages


if(length(ip)) devtools::install_github(rp2, dependencies = TRUE) # install packages


lapply(rp, require, character.only = TRUE) # load packages

rm(rp, ip, rp2) # remove vector
### Directory setup ###########################################################

setwd("D:/EAGLE/Programming/PRoject")

outputfolder <- "D:/EAGLE/Programming/PRoject"

# optional if more GDAL versions installed
gdal_setInstallation(search_path = "C:\\OSGeo4W64\\bin")



# SAGA setup (required for terrain analysis)

# path to SAGA folder with saga.exe
path <- "D:\\EAGLE\\Geoanalysis\\SAGA\\saga-8.1.1_x64\\saga-8.1.1_x64"

saga <- saga_gis(saga_bin = paste0(path,"\\saga_cmd.exe"), raster_format = "GeoTIFF")


### define species and AOI #####################################################

### latin scientific species name
species_name <- "arnica montana"

# processing for entire country takes a while, here selection of the Austrian Bundesland Vorarlberg
country <- "Austria"
level <- 1

shapeAOI <- getData(name = "GADM", country= country, level =level)
shapeAOI@data[["NAME_1"]]
shapeAOI <- shapeAOI[shapeAOI$NAME_1 == "Vorarlberg",]
shapeAOI <- st_as_sf(shapeAOI[, -(1:ncol(shapeAOI))])

plot(shapeAOI)

### load your own shapefile 
#shapeAOI <- st_read(path.shp")


### Load data #################################################################
### GADM adminstrative data




### reproject to UTM and setting up variables

  if (!st_crs(shapeAOI)== 4326){
    shapeAOI <- st_transform(shapeAOI,4326)
  } 
shapeAOI_wkt <- sf_convert(shapeAOI)
bbox_ll <- extent(shapeAOI)
coords <- st_centroid(shapeAOI) %>% 
                st_coordinates()
UTM_zone <- floor((coords[1,] + 180) / 6) + 1
print(paste0("UTM zone is: ", UTM_zone[1]))
crs_UTM  <-  paste0( "+proj=utm +zone=", UTM_zone[1], " +datum=WGS84 +units=m +no_defs")
crs_WGS84 <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" 
shapeAOI_UTM <- st_transform(shapeAOI, crs =crs_UTM)
   
crs_igh='+proj=igh +lat_0=0 +lon_0=0 +datum=WGS84 +units=m +no_defs'
extent_AOI_UTM <- extent(st_buffer(shapeAOI_UTM[,1:(ncol(shapeAOI_UTM)-1)], dist = 2000))
extent_AOI_igh <- extent(st_transform(st_buffer(shapeAOI_UTM[,1:(ncol(shapeAOI_UTM)-1)], dist = 2000), crs=crs_igh))





### GBIF species data #########################################################
# https://www.r-bloggers.com/2021/03/downloading-and-cleaning-gbif-data-with-r/


gbif_data <- occ_data(scientificName = species_name, hasCoordinate = TRUE, geometry = shapeAOI_wkt, geom_big = "bbox", limit = 100000)
data_citation <- gbif_citation(gbif_data)
print(data_citation)

names(gbif_data$data)
species_coords <- gbif_data$data[ , c("decimalLongitude", "decimalLatitude", "individualCount", "occurrenceStatus", "coordinateUncertaintyInMeters", "institutionCode", "references", "year", "month", "day")] %>% date_create(year, month, day) 
species_coords <- species_coords[,!names(species_coords) %in% c("year","month", "day")]
head(species_coords)

# cleaning data
nrow(species_coords)
species_coords <- species_coords[!(species_coords$occurrenceStatus %in% c("absent", "Absent", "ABSENT", "ausente", "Ausente", "AUSENTE")),]

vague_rows <- which(!is.na(species_coords$coordinateUncertaintyInMeters) & species_coords$coordinateUncertaintyInMeters > 1000)
species_coords <- species_coords[-vague_rows,]
nrow(species_coords)

species_coords <- coord_incomplete(coord_imprecise(coord_impossible(coord_unlikely(species_coords))))
nrow(species_coords)

# duplicates removal, sampling Bias reduction
species_coords <- dedup(species_coords)
nrow(species_coords)

species_points <- st_as_sf(species_coords, coords = c("decimalLongitude", "decimalLatitude"), crs=crs_WGS84) %>% st_transform(crs_UTM)

plot(shapeAOI_UTM)
plot(species_points, add = T)





### ISRIC SOIL Data ###########################################################

# https://www.isric.org/explore/soilgrids
# https://git.wur.nl/isric/soilgrids/soilgrids.notebooks/-/blob/master/markdown/wcs_from_R.md
# https://git.wur.nl/isric/soilgrids/soilgrids.notebooks/-/blob/master/markdown/webdav_from_R.md
# https://rpubs.com/ials2un/soilgrids_webdav # that is the one implemented here


dir.create("soildata")

### boundiing box in interrupted Goode homolosine projection, it is better handled by the webservice

bb=c(round(as.numeric(extent_AOI_igh@xmin), 0), round(as.numeric(extent_AOI_igh@ymax), 0),round(as.numeric(extent_AOI_igh@xmax), 0),round(as.numeric(extent_AOI_igh@ymin), 0)) 

soilgriddata <- c("bdod", "cec", "cfvo", "clay", "nitrogen", "phh2o", "sand", "silt", "soc", "ocd")



for (i in 1:length(soilgriddata)){
  sg_url <- "/vsicurl/https://files.isric.org/soilgrids/latest/data/"
  soildata  <-  paste0(soilgriddata[i], "/", soilgriddata[i], "_15-30cm_mean.vrt")
  soiloutput <- paste0("./soildata/", soilgriddata[i], "_15_30.tif")
  gdal_translate(paste0(sg_url,soildata), soiloutput ,
                 tr=c(250,250),
                 projwin=bb,
                 projwin_srs =crs_igh,
                 verbose=TRUE)
  soildata2  <-  paste0(soilgriddata[i], "/", soilgriddata[i], "_5-15cm_mean.vrt")
  soiloutput2 <- paste0("./soildata/", soilgriddata[i], "_5_15.tif")
  gdal_translate(paste0(sg_url,soildata2), soiloutput2 ,
                 tr=c(250,250),
                 projwin=bb,
                 projwin_srs =crs_igh,
                 verbose=TRUE)
  soildata3  <-  paste0(soilgriddata[i], "/", soilgriddata[i], "_0-5cm_mean.vrt")
  soiloutput3 <- paste0("./soildata/", soilgriddata[i], "_0_5.tif")
  gdal_translate(paste0(sg_url,soildata3), soiloutput3 ,
                 tr=c(250,250),
                 projwin=bb,
                 projwin_srs =crs_igh,
                 verbose=TRUE)
  rm(soildata, soiloutput,soildata2, soiloutput2,soildata3, soiloutput3, sg_url, i)
}


soil_files <- list.files("./soildata/", full.names = T, pattern = "*.tif")

dir.create("soildata/warped")

for (i in 1:length(soil_files)){
    if (i == 1){
    gdalwarp(soil_files[i], paste0("./soildata/warped", str_sub(soil_files[i], start= 11, end = -5), "_warped.tif"), s_srs =crs_igh, t_srs = crs_UTM, tr= c(xres=250, yres=250), overwrite = T)
    soil_raster <- raster(paste0("./soildata/warped", str_sub(soil_files[i], start= 11, end = -5), "_warped.tif"))
  } 
  else {
    gdalwarp(soil_files[i], paste0("./soildata/warped", str_sub(soil_files[i], start= 11, end = -5), "_warped.tif"), s_srs =crs_igh, t_srs = crs_UTM, tr= c(xres=250, yres=250), overwrite = T)
    soil_raster <- stack(soil_raster, raster(paste0("./soildata/warped", str_sub(soil_files[i], start= 11, end = -5), "_warped.tif")))
  }
}

soil_names <- c("bulk_density_0_5", "bulk_density_15_30", "bulk_density_5_15", "cec_0_5" ,  "cec_15_30",  "cec_5_15",  "fraction_of_coarse_fragments_0_5", "fraction_of_coarse_fragments_15_30", "fraction_of_coarse_fragments_5_15", "clay_0_5", "clay_15_30", "clay_5_15", "nitrogen_0_5", "nitrogen_15_30", "nitrogen_5_15", "organic_carbon_density_0_5",  "organic_carbon_density_15_30",  "organic_carbon_density_5_15",  "phh2o_0_5", "phh2o_15_30", "phh2o_5_15", "sand_0_5", "sand_15_30", "sand_5_15", "silt_0_5", "silt_15_30", "silt_5_15", "soil_organic_content_0_5", "soil_organic_content_15_30", "soil_organic_content_5_15")
names(soil_raster) <- soil_names
res(soil_raster)

plot(soil_raster$bulk_density_0_5)
plot(shapeAOI_UTM, add =T)
plot(species_points, add = T)



### soilgrids stores data in integer format therefore conversion in conventional units

for (i in 1:nlayers(soil_raster)) {
  if (names(soil_raster[[i]]) %in% grep("bulk_density*", soil_names, value =T) | names(soil_raster[[i]]) %in% grep("nitrogen*", soil_names, value =T) ) {
 soil_raster[[i]] <- soil_raster[[i]]/100
} else {
  soil_raster[[i]] <- soil_raster[[i]]/10
}
}




### CORINE Landcover ###########################################################
### extent only Europe

Corine <- raster(paste0(outputfolder, "/U2018_CLC2018_V2020_20u1.tif")) 
Corine <- projectRaster(Corine, soil_raster, method = "ngb" )
names(Corine) <- "landcover2018"


### explanation to classes, apart from that csv incorporated
clc_classes <- read.csv(paste0(outputfolder, "/clc_legend.csv"), header = T)
print(clc_classes)

### categorical raster data
Corine <- as.factor(reclassify(Corine, matrix(c(44, Inf, NA), ncol = 3, byrow = T)))

### masks for displaying
water_mask <- rasterToPolygons(mask(reclassify(Corine, matrix(c(0, 39, NA,  39, 44, 1,  45, Inf, NA), ncol = 3, byrow=T)), shapeAOI_UTM), dissolve = T)
urban_mask <- rasterToPolygons(mask(reclassify(Corine, matrix(c(1, 6, 1,  6, 45, NA), ncol = 3, byrow=T)), shapeAOI_UTM), dissolve = T)


plot(water_mask)
plot(urban_mask)



### DEM #######################################################################

### define download SRTM extent (might not work for large countries)
### trying to ensure that entire AOI is downloaded as SRTM tiles using the corner points and centroid, alternatively adaption working with SRTM Tiles .shp feasible

SRTM_Download_for_cornerpoints <- function (boundingbox, centroid_coords) {
  SW <- st_point(c(boundingbox@xmin, boundingbox@ymin))
  SE <- st_point(c(boundingbox@xmax, boundingbox@ymin))
  NW <- st_point(c(boundingbox@xmin, boundingbox@ymax))
  NE <- st_point(c(boundingbox@xmax, boundingbox@ymax))
  SRTM <- getData("SRTM", lon = centroid_coords[1], lat = centroid_coords[2] )
  SRTM_extent <- st_as_sf(as(extent(SRTM), "SpatialPolygons"))
  # check if any corner points are not contained and if so download another tile  
    
  if (!st_contains(SW,SRTM_extent, sparse = F)) { 
  SRTM2 <- getData("SRTM", lon = boundingbox@xmin, lat = boundingbox@ymin)
  SRTM_list <- list(SRTM, SRTM2)
  SRTM_list$fun <- mean
  SRTM <- do.call(mosaic, SRTM_list)
  SRTM_extent <- st_as_sf(as(extent(SRTM), "SpatialPolygons"))
  rm(SRTM_list, SRTM2)
  }
  if (!st_contains(SE,SRTM_extent, sparse = F)) { 
    SRTM2 <- getData("SRTM", lon = boundingbox@xmax, lat = boundingbox@ymin)
    SRTM_list <- list(SRTM, SRTM2)
    SRTM_list$fun <- mean
    SRTM <- do.call(mosaic, SRTM_list)
    SRTM_extent <- st_as_sf(as(extent(SRTM), "SpatialPolygons"))
    rm(SRTM_list, SRTM2)
  } 
  if (!st_contains(NW, SRTM_extent, sparse = F)) { 
    SRTM2 <- getData("SRTM", lon = boundingbox@xmin, lat = boundingbox@ymax)
    SRTM_list <- list(SRTM, SRTM2)
    SRTM_list$fun <- mean
    SRTM <- do.call(mosaic, SRTM_list)
    SRTM_extent <- st_as_sf(as(extent(SRTM), "SpatialPolygons"))
    rm(SRTM_list, SRTM2)
  }
  if (!st_contains(NE, SRTM_extent, sparse = F)) { 
    SRTM2 <- getData("SRTM", lon = boundingbox@xmax, lat = boundingbox@ymax)
    SRTM_list <- list(SRTM, SRTM2)
    SRTM_list$fun <- mean
    SRTM <- do.call(mosaic, SRTM_list)
    rm(SRTM_list, SRTM2, SRTM_extent)
  }
  return(SRTM)
}

SRTM <- SRTM_Download_for_cornerpoints(bbox_ll, coords)

plot(SRTM)
plot(shapeAOI, add=T)

dir.create("terrain")
SRTM <- projectRaster(SRTM, soil_raster, method = "bilinear", filename = "./terrain/DEM.tif", overwrite=T ) 
plot(SRTM) 




### calculation of terrain parameters with SAGA ###############################

setwd("./terrain")

saga$ta_compound$basic_terrain_analysis(elevation = SRTM,
                                        slope = "slope.tif",
                                        aspect = "aspect.tif",
                                        hcurv = "plan_curv.tif",
                                        vcurv = "profile_curv.tif",
                                        convergence = "convergence.tif",
                                        flow = "catchment_area.tif",
                                        wetness = "TWI.tif",
                                        lsfactor = "LSfactor.tif",
                                        chnl_base = "channel_base.tif",
                                        chnl_dist = "channel_distance.tif",
                                        vall_depth = "valley_depth.tif",
                                        rsp = "relative_slope_position.tif")




saga_raster_files <- list.files(paste0(outputfolder,"\\terrain"), full.names = T,  pattern="*.tif$") 


### building stack of rasters
for (i in 1:length(saga_raster_files)){
  if (i == 1){
    terrain_raster <- raster(saga_raster_files[i])
  }
  else {
    terrain_raster <- stack(terrain_raster,raster(saga_raster_files[i]))
    }
}  

print(terrain_raster) 


#adjust slope and Aspect 
terrain_raster$slope <- (terrain_raster$slope)*180/pi 
terrain_raster$aspect <- (terrain_raster$aspect)*180/pi 
terrain_raster$relative_slope_position <- (terrain_raster$relative_slope_position)*180/pi 

terrain_raster$eastness <- sin(terrain_raster$aspect)
terrain_raster$northness <- cos(terrain_raster$aspect)

terrain_raster <- dropLayer(terrain_raster, match("aspect", names(terrain_raster)))
setwd("..")





### Tree Biomass data ##########################################################
### Download data from CEDA does not work automatically yet. Required is a certificate that can be obtained via a python script. 
### However incorporation in R script did not work out yet. 
### Data available at: https://data.ceda.ac.uk/neodc/esacci/biomass/data/agb/maps/v2.0/geotiff/2018
### citation: Santoro, M.; Cartus, O. (2021): ESA Biomass Climate Change Initiative (Biomass_cci): Global datasets of forest above-ground biomass for the years 2010, 2017 and 2018, v3. NERC EDS Centre for Environmental Data Analysis, 26 November 2021. doi:10.5285/5f331c418e9f4935b8eb1b836f8a91b8


biomass_raw <- list.files(paste0( "./biomass"), full.names = T,  pattern="*.tif$", recursive=T) 

biomass.list <- lapply(1:length(biomass_raw), function(x) {raster(biomass_raw[x]) })

biomass.list$fun <- mean
biomass.mosaic <- do.call(mosaic,biomass.list) %>% projectRaster(soil_raster, method = "bilinear")
names(biomass.mosaic) <- "biomass"
res(biomass.mosaic)
plot(biomass.mosaic)





### World Clim Data ###########################################################


Worldclim <- getData("worldclim", var = "bio", res=0.5, lon=round(coords[1]), lat= round(coords[2])) %>% 
  projectRaster(soil_raster, method = "bilinear")

World_clim_names <- c("BIO1_mean_temp", "BIO2_meanDiurnal_temp", "BIO3_iso_Temp", "BIO4_Seasonality_Temp", "BIO5_max_temp", "BIO6_min_temp", "BIO7_range_temp", "BIO8_mean_temp_wet", "BIO9_mean_temp_dry", "BIO10_mean_temp_warm", "BIO11_mean_temp_cold", "BIO12_annPrec", "BIO13_prec_wetmonth", "BIO14_prec_drymonth", "BIO15_prec_seasonality", "BIO16_prec_wetquart", "BIO17_prec_dryquart", "BIO18_prec_warmquart", "BIO19_prec_coldquart")
names(Worldclim) <- World_clim_names

Worldclim

plot(Worldclim$BIO1_mean_temp)
plot(shapeAOI_UTM, add=T)

# conversion of temperature (due to online integer storage convention)

for (i in 1:nlayers(Worldclim)) {
  if (names(Worldclim[[i]]) %in% grep("temp*", World_clim_names, value =T)) {
    Worldclim[[i]] <- Worldclim[[i]]/10
  } else{}
  if (names(Worldclim[[i]]) %in% grep("Temp*", World_clim_names, value =T)) {
   Worldclim[[i]] <- Worldclim[[i]]/100
  }
}






### extract covariates data to points ##########################################
# extraction for corine landcover as factor (without interpolation), bilinear for metric scales

covariates <- stack(soil_raster, terrain_raster, biomass.mosaic, Worldclim, Corine)


species_present <- data.frame(species_points$geometry, raster::extract(covariates[[1:(nlayers(covariates)-1)]], species_points, method = "bilinear"), raster::extract(covariates$landcover2018, species_points, method = "simple"))

# name change of landcover column during extraction
colnames(species_present)[ncol(species_present)] <- "landcover2018"


nrow(species_present)
species_present_omit <- na.omit(species_present)
nrow(species_present_omit)

# species occurences that fall on NA pixels (for final plot)
species_presence_na <- st_as_sf(species_present[rowSums(is.na(species_present)) > 0, ], crs=crs_UTM) 
plot(species_presence_na$geometry)


### create background data #####################################################
# non NA values raster mask to create points only in non NA area

covariates_mask <- !is.na(sum(mask(covariates, shapeAOI_UTM)))
covariates_mask <- reclassify(covariates_mask, c(-Inf, 0, NA), ncol = 3, byrow=T)

plot(covariates_mask)


# based on non NA mask random points creation
# here by polygonizing raster, quite cumbersome
poly_raster <- rasterToPolygons(covariates_mask, fun = function (x) {x==1}, na.rm=T, dissolve = T)
poly_raster
set.seed(0)
background_points <- st_as_sf(spsample(poly_raster[poly_raster$layer ==1], nrow(species_present_omit),  type="random"), crs=crs_UTM)

plot(poly_raster, add=T)
plot(background_points, col = "orange", add = T )

# alternative with contour line polygonize faster but doesn't work yet
# see: https://newbedev.com/generating-random-points-inside-raster-boundary-using-r


### calculate point density per raster cells
point_density <- nrow(species_present_omit)/cellStats(covariates_mask, 'sum')
point_density

### extract raster values for background data and create complete dataframe
background <- data.frame(background_points$geometry, raster::extract(covariates[[1:(nlayers(covariates)-1)]], background_points, method = "bilinear"), raster::extract(covariates$landcover2018, background_points, method = "simple"))

# name change of landcover column during extraction
colnames(background)[ncol(background)] <- "landcover2018"

pres <- as.factor(c(rep(1, nrow(species_present_omit)), rep(0, nrow(background))))
sdm_data <- data.frame(cbind(pres, rbind(species_present_omit, background)))

head(sdm_data)

### multicollinearity ##########################################################


### correlation between terrain_rasters

corrplot_calculation <- function (raster_stack) {
  cor <- layerStats(raster_stack, 'pearson', na.rm=T)
  cor <- cor[[1]]
  cor[which(cor>1)] <- 1
  cor[which(cor<(-1))] <- -1
  stackcor <- corrplot(cor, method = "circle", type = "lower", order = "FPC", col = rev.default(COL2('RdBu', 200)), tl.col = "#010101", tl.srt = 45)
  return(stackcor)
}

terrain_corplot <- corrplot_calculation(terrain_raster)
soil_corplot <- corrplot_calculation(soil_raster)
climate_corplot <- corrplot_calculation(Worldclim)


### testing for multicollinearity 


mcl <- multi.collinear(sdm_data[,3:ncol(sdm_data)], p=0.05, na.rm = T)
print(mcl)

sdm_data <- sdm_data[,-which(names(sdm_data) %in% mcl )]






### Random Forest Modelling ####################################################

### mostly done according to Jeffrey Evans Tutorial: https://evansmurphy.wixsite.com/evansspatial/random-forest-sdm

# create random forest model parameters
rf_model <- rf.modelSel(x=sdm_data[,3:ncol(sdm_data)], y=sdm_data[,"pres"], imp.scale="mir", ntree=1500, seed = 12)
rf_model$selvars

# select the important covariates

sel_vars <- rf_model$selvars

# Train the model
rf_fit <- randomForest(y=sdm_data[, "pres"], x=sdm_data[,sel_vars], ntree=1500, importance=TRUE, norm.votes=TRUE, proximity=TRUE)
rf_fit
plot(rf_fit, main="Bootstrap Error Convergence")

rf_covariates <- covariates[[rownames(rf_fit$importance)]]

# make prediction on entire AOI
rf_prediction <- predict(rf_covariates, rf_fit, filename= paste0(outputfolder, "/", species_name, "_RF_model.grd"), type="prob", index=2, na.rm=TRUE, overwrite=TRUE, progress="window")
rf_prediction


palette = brewer.pal(9,'YlGn')
palette_water = "dodgerblue"


pdf(paste0("SDM_", species_name), width = 14, height = 10)
plot(rf_prediction, col=palette, main= species_name)
plot(shapeAOI_UTM, add=T)
plot(water_mask, col= "blue", add = T)
plot(urban_mask, col= "black", add = T)
plot(species_points, add=T, col="red", pch=20, cex=0.75)
plot(background_points, add=T, col="purple", pch=20, cex=0.75)
plot(species_presence_na, add=T, pch= 24, cex=1.2, col= "yellow")
dev.off()


### covariates importance
pdf("covariates_importance", width = 8, height= 10)
RF_importance_matrix <- as.matrix(rf_fit$importance[,3])   
ord <- rev(order(RF_importance_matrix[,1], decreasing=TRUE)[1:dim(RF_importance_matrix)[1]]) 
dotchart(RF_importance_matrix[ord,1], main="Scaled Variable Importance", cex=0.8, pch=19) 
dev.off()

### partial plot

pdf("partial_plot", width = 12, height = 12)
par(mfrow=c(2,2))
ord2 <- rev(ord)
for(i in ord2[1:4]) {
 rf.partial.prob(rf_fit, sdm_data[,sel_vars], i, "1", smooth="spline", xlab=sel_vars[i], raw.line=FALSE)
}  
dev.off()

### dimension matrix
rf_cmd <- cmdscale(1 - rf_fit$proximity, eig=TRUE, k=5) 
pa_col=ifelse(sdm_data[,"pres"] == "1", "red", ifelse(sdm_data[,"pres"] == "0", "black", NA))
rf_cmd <- data.frame(rf_cmd$points)

pdf("dimension_matrix", width = 10, height = 10)
plot(rf_cmd[,1:2], ylab="DIM 1", xlab="DIM 2", pch=16, col=pa_col, main= paste("Pres/Abs", "PROXIMITY MATRIX MDS d=2", sep=" - "))
legend("topright",pch=c(16,16), col=c("black", "red"),
       legend=c("Absent","Present"))
dev.off()





### Model validation ###########################################################
### Bootstrapping
rf_perm <- rf.significance(rf_fit, sdm_data[,sel_vars], nperm = 399, ntree = 1501)
rf_perm
plot(rf.perm)

### Cross Validation
rf_cv <- rf.crossValidation(rf_fit, sdm_data[,sel_vars], seed=12, p=0.10, n=399, ntree=1501)
print(rf_cv)

### future tasks: include different models, e.g. utilizing biomod2 or dismo package
                  # more covariates for automatic download

