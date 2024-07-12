# 02_Variables.R
# Step 2

## In this step, we must download all relevant rasters, and crop them for uses in later scripts

# Tori Ford
# March 27 2024

library(geodata)
library(dplyr)
library(raster)


# Download ----------------------------------------------------------------

# Create folder for rasters
dir.create("raw/rast")


# World Clim: https://www.worldclim.org/data/bioclim.html

# ## Create a list for sets of variables downloaded with same naming convention, update as you go
wc_var <- c("bio","srad","vapr")

## Download for WorldClim Data

# Open Empty Lists, this is the option if you don't want to download, see line 34
clim <- c()

## Loop for WorldClim Data by cyclng variables
for (q in 1:length(wc_var)) {
	var <- wc_var[q]
	print(var)
	
	worldclim <- geodata::worldclim_global(var = var,
																				 res = 2.5, # Select an appropriate resolution, this matters if it is not 30s, see line 42
																				 path = paste0("raw/rast/")
																				 # path = tempdir() # Files are temporarily stored in R
																				 )
	# clim[[q]] <- worldclim
	}

# IF you choose a resolution besides 30s, you will have to resample the soil rasters to match the resolution between predictors
# in this case, we are working in 2.5 arc-minutes, soil rasters are 30s so they must be resampled.
# Read in a bioclim variable as reference
one.samp <- raster("raw/rast/XX.tif")

# Use to STOP aux files when using writeRaster function
setGDALconfig("GDAL_PAM_ENABLED", "FALSE") 


# Load World Soils Data


# List environmental soil predictors of interest
sworld_var <- c("sand", "silt", "soc", "bdod","cfvo", "clay", "nitrogen","ocd","phh2o")  

# Set up variables + depths: check available depths/variables in the geodata vignette: https://cran.r-project.org/web/packages/geodata/geodata.pdf
five <- paste(sworld_var,"_5",sep = "")
fifteen <- paste(sworld_var,"_15",sep = "")
thirty <- paste(sworld_var,"_30",sep = "")

# Total list of variables + depths
soils <- c(five,fifteen,thirty)

# Open Empty Lists, this is the option if you don't want to download, see line 74
# soil <- c()

for (p in 1:length(soils)) {

	s <- strsplit(soils[[p]],"[_]")[[1]][1] # Split variable name for soil_world(var)
	d <- as.numeric(strsplit(soils[[p]],"[_]")[[1]][2]) # Split variable name for soil_world(depth)
	print(paste0(s,d))

	soilworld <- geodata::soil_world(
		var = s,
		depth = d,
		stat="mean",
		path = tempdir() # Files are temporarily stored in R
		# path = "raw/rast" # option to download the 30s rasters
	)
	soil.res <- terra::resample(soilworld, one.samp)
	# soil[[p]] <- soil.res # IF you don't wnt to download the resampled rasters, open list on 58, and use this line
	writeRaster(soil.res, paste0("raw/rast/",soils[p], ".grd")) # option to download resampled rasters

	
}


## IF DOWNLOADED, read in the rasters
soil.ras <- stack(list.files(paste0("raw/rast/"), pattern = "*.tif", full.names = TRUE))
bioclim.ras <- stack(list.files(paste0("raw/rast/"), pattern = "*.tif", full.names = TRUE))

## Extend Rasters to Match Smaller Extent (soil, typically)
beginCluster(6)

## Open empty list for cropped variabels
bioclim.crop <- c()

for (c in 1:length(names(bioclim.ras))) {
	
	# Subset layer
	layer <- bioclim.ras[[c]]
	
	clim.crop <- terra::crop(layer,ext(soil.ras))
	
	bioclim.crop[[c]] <- clim.crop
	
}
endCluster()

# Combine them as predictors
predictors <- c(soil.ras, bioclim.crop)

# Consider saving the global rasters at this point


# ## IF SAVED TEMP, unlist the rasters

# soil.ras <- rast(soil)
# bioclim.ras <- rast(clim)
# 
# ## Extend Rasters to Match Smaller Extent (soil, typically)
# beginCluster(6)
# 
# ## Open empty list for cropped variabels
# bioclim.crop <- c()

# for (c in 1:length(names(bioclim.ras))) {
# 	
# 	# Subset layer
# 	layer <- bioclim.ras[[c]]
# 	
# 	clim.crop <- terra::crop(layer,ext(soil.ras))
# 	
# 	bioclim.crop[[c]] <- clim.crop
# 	
# }
# endCluster()
# 
# # Combine them as predictors
# predictors <- c(soil.ras, bioclim.crop)

# Crop to Total (study) Extent ----------------------------------------------------


## Read data in; occurrences will be used to constrain occurrences
datum <- read.csv("raw/CJG_popK2_for_TF.csv")

## Change colum names
datum <- datum %>% rename(
	lon = LON,
	lat = LAT
)

## Set up buffer for crop + mask
## Total Extent
total.sf <- sf::st_as_sf(datum, coords = c("lon","lat"), crs = terra::crs(predictors))


## Get bbox for raster crop to total study extent (ALL POPs)
total.ext <- total.sf %>% 
	st_as_sf(coords = c("x","y"), crs = 4326) %>% 
	st_bbox()



## Crop all predictors by outer bounding box limits
bbox.pred <- raster::crop(predictors, extent(ssp.ext)) # crop by extent
bbox.pred <- raster::mask(bbox.pred, extent(ssp.ext)) # mask cropped values with NAs
for (b in 1:length(names(bbox.pred))) {
	layer <- bbox.pred[[b]] # Subset out layer
	name <- names(layer) # Get layer name
	raster::writeRaster(layer, paste0("raw/rast/total/",name,".asc"), overwrite = T) #write layer out with matching name
	print(Sys.time()) # print time taken per write for records, comment out if not needed. 
}


# Crop to Each Populations Accessible Area (Buffers) ----------------------

## Read data in; occurrences will be used to constrain occurrences
ssp <- read.csv("raw/CJG_popK2_for_TF.csv")

## Change colum names
spp <- spp %>% rename(
	lon = LON,
	lat = LAT
)

# For "buffer"=b  in 1:n(unique populations in the dataset)
for (b in unique(spp$POPk2)) {
	
	## Species Extent for ENM
	## THIS IS THE BUFFERS FOR THE MODEL ITSELF
	## Split Occurences Per K Group, One total df, one coordinates
	oocs <- ssp[grep(b, ssp$POPk2),]
	xy <- oocs[, c("LON", "LAT")]
	xy.sf <- sf::st_as_sf(xy, coords = c("LON","LAT"), crs = raster::crs(bbox_envs))
	# xy.sf <- sf::st_as_sf(xy, coords = c("LON","LAT"), crs = raster::crs(predictors))
	
	
	## Convert projection to equal area for translation to meters
	eckertIV <- "+proj=eck4 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"
	xy.sf <- sf::st_transform(xy.sf, crs = eckertIV)
	
	# Buffer all occurrences by 500 km, union the polygons together 
	# (for visualization), and convert back to a form that the raster package 
	# can use. Finally, we reproject the buffers back to WGS84 (lat/lon).
	# We choose 100 km here per: https://www.scielo.br/j/mioc/a/d45vRjYQKv96GDRGrWLsh9x/?lang=en&format=html 
	xy.buf <- sf::st_buffer(xy.sf, dist = 100000) %>% 
		sf::st_union() %>% 
		sf::st_sf() %>%
		sf::st_transform(crs = raster::crs(bbox_envs))
}
## Species Extent for ENM
## THIS IS THE BUFFERS FOR THE MODEL ITSELF
## Split Occurences Per K Group, One total df, one coordinates
oocs <- ssp[grep(s, ssp$POPk2),]
xy <- oocs[, c("x", "y")]
xy.sf <- sf::st_as_sf(xy, coords = c("x","y"), crs = raster::crs(envs))


## Convert projection to equal area for translation to meters
eckertIV <- "+proj=eck4 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"
xy.sf <- sf::st_transform(xy.sf, crs = eckertIV)

# Buffer all occurrences by 500 km, union the polygons together 
# (for visualization), and convert back to a form that the raster package 
# can use. Finally, we reproject the buffers back to WGS84 (lat/lon).
# We choose 100 km here per: https://www.scielo.br/j/mioc/a/d45vRjYQKv96GDRGrWLsh9x/?lang=en&format=html 
xy.buf <- sf::st_buffer(xy.sf, dist = 100000) %>% 
	sf::st_union() %>% 
	sf::st_sf() %>%
	sf::st_transform(crs = raster::crs(envs))



dir.create(paste0("raw/rast/buf",b))

## CROP, MASK, SAVE RADIAL BUFFER
for (j in 1:length(rast_files)) {
	
	# Subset raster layer
	layer <- predictors[[j]]
	
	# Setup file names
	name <- names(layer)
	
	print(name) ## Just here to check on progress of loop
	
	# Crop and mask
	cropped <- terra::crop(layer, xy.buf) #ext
	masked <- terra::mask(cropped, xy.buf)
	
	# Write raster
	terra::writeRaster(masked,
										 filename = paste0("raw/rast/buf",b,"/",name,".asc"),
										 overwrite = F)
	print(paste0(j))
}
