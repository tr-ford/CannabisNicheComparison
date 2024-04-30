# Variables.R
# Step 2

# Tori Ford
# March 27 2024

library(geodata)
library(dplyr)
library(raster)


# Download ----------------------------------------------------------------


# Download World Soils Data

sworld_var <- c("sand", "silt", "soc", "bdod","cfvo", "clay", "nitrogen","ocd","phh2o")
depth <- c(5, 15, 30) # , 100, 200 cm not working atm, prolly use diff script

## Create SoilGridsV2 World Soils Directory
# dir.create("raw/rasters/worldsoils")


for (p in 1:length(sworld_var)) {
	for (d in depth) {

		print(d)
		name <- sworld_var[p]
		print(name)

		soilworld <- geodata::soil_world(
			var = name,
			depth = d,
			stat="mean",
			path = "raw/rast"
		) ## Soil download not working, use VSI to interface in the meantime, read to local variable to resample
	}
}


# ## Create a list for sets of variables downloaded with same naming convention, update as you go
wc_var <- c("bio","srad","vapr")

## Download for WorldClim Data

## Loop for WorldClim Data
for (q in 1:length(wc_var)) {
	var <- wc_var[q]
	print(var)
	
	worldclim <- geodata::worldclim_global(var = var,
																				 res = 2.5,
																				 path = paste0("raw/rast/",var))
	
}


# Crop to Study Extent ----------------------------------------------------


## Read data in; occurrences
datum <- read.csv("raw/CJG_popK2_for_TF.csv")

## Change colum names
datum <- datum %>% rename(
	lon = LON,
	lat = LAT
)

# hist_clim <- worldclim_global(var = "bio", res = 10, version="2.1", path=tempdir())

## List clim files, read them in together
list1 <- list.files(paste0("raw/rast/bio/"), pattern = "*.tif", full.names = TRUE) 
# list4 <- list.files(paste0("raw/rast/soil_world/"), pattern = "*.asc", full.names = TRUE) 
list2 <- list.files("raw/rast/srad/", pattern = "*.tif", full.names = TRUE) 
list3 <- list.files("raw/rast/vapr/", pattern = "*.tif", full.names = TRUE)

lists <- c(list1,list2,list3)

## Set up buffer for crop + mask

# We'll now experiment with a different spatial R package called sf (simple features).
# Let's make our occs into a sf object -- as the coordinate reference system (crs) for these 
# points is WGS84, a geographic crs (lat/lon) and the same as our envs rasters, we specify it 
# as the RasterStack's crs.
datum.sf <- sf::st_as_sf(datum, coords = c("lon","lat"), crs = raster::crs(rsclat))

# Now, we project our point data to an equal-area projection, which converts our 
# degrees to meters, which is ideal for buffering (the next step). 
# We use the typical Eckert IV projection.
eckertIV <- "+proj=eck4 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"
datum.sf <- sf::st_transform(datum.sf, crs = eckertIV)

# Buffer all occurrences by 500 km, union the polygons together 
# (for visualization), and convert back to a form that the raster package 
# can use. Finally, we reproject the buffers back to WGS84 (lat/lon).
# We choose 500 km here to avoid sampling the Caribbean islands.
datum.buf <- sf::st_buffer(datum.sf, dist = 50000) %>% 
	sf::st_union() %>% 
	sf::st_sf() %>%
	sf::st_transform(crs = raster::crs(rsclat))
plot(rsclat, main = names(rsclat)[1])
points(datum)
# To add sf objects to a plot, use add = TRUE
plot(datum.buf, border = "blue", lwd = 3, add = TRUE)



## Crop to extents of shape file
for (o in 1:length(lists)) {
	# Subset raster layer
	layer <- rast(lists[[o]])
	# Setup file names
	name <- names(layer)
	print(name) ## Just here to check on progress of loop
	# out <- paste0(path, name)
	# outfile <- paste0(out, end)
	# Crop and mask
	cropped <- terra::crop(layer, datum.buf) #ext
	masked <- terra::mask(cropped, datum.buf)
	# Write raster
	terra::writeRaster(masked,
										 filename = paste0("raw/rast/buf/",name,".asc"),
										 overwrite = F)
	print(paste0(o))
}


lists <- list.files(paste0("raw/rast/buf/"), pattern = "*.grd", full.names = TRUE) 
## Converts buffered rasters from .grd to .asc for MAXENT 
for (a in 1:length(lists)) {

	laya <- trim(rast(lists[[a]]))
	laya[is.nan(laya)] <- 0
	laya[is.na(laya)] <- 0
	name <- names(laya)
	terra::writeRaster(laya, paste0("raw/rast/buf/",name,".asc"), overwrite = T)
}
