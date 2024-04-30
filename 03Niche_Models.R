# Niche_Models.R
# Step 4

# Tori Ford
# April 2, 2024
# Objective: Complete All Steps of Ecological NIche Models, From Maxent
# https://rpubs.com/janikh99/1044769
# Modified from tutorial, native dynamics vs invaded



library(ecospat)
library(geodata)
library(raster)
# library(maptools)
library(sf)
library(dplyr)
library(devtools)
library(terra)
library(ade4)
library(ape)
library(biomod2)
library(spam)
library(spam64)


# Premodel ----------------------------------------------------------------


datum <- read.csv("raw/CJG_popk2_for_TF.csv")
datum <- datum[,c(3,4)]
s_datum <- read.csv("raw/Senv1.csv")
s_datum <- s_datum[,c(1,2)]
s_datum <- s_datum[complete.cases(data.frame(s_datum)), ]
n_datum <- read.csv("raw/Nenv2.csv")
n_datum <- n_datum[,c(1,2)]
n_datum <- n_datum[complete.cases(data.frame(n_datum)), ]

# # Load World Soils Data
# 
# sworld_var <- c("sand", "silt", "soc", "bdod","cfvo", "clay", "nitrogen","ocd","phh2o")
# depth <- c(5, 15, 30) 
# 
# # Set up variables
# five <- paste(sworld_var,"_5",sep = "")
# fifteen <- paste(sworld_var,"_15",sep = "")
# thirty <- paste(sworld_var,"_30",sep = "")
# 
# # Total list of variables + depths
# soils <- c(five,fifteen,thirty)
# 
# # Open Empty Lists
# soil <- c()
# 
# for (p in 1:length(soils)) {
# 	
# 		s <- strsplit(soils[[p]],"[_]")[[1]][1]
# 		d <- as.numeric(strsplit(soils[[p]],"[_]")[[1]][2])
# 		print(paste0(s,d))
# 		
# 		soilworld <- geodata::soil_world(
# 			var = s,
# 			depth = d,
# 			stat="mean",
# 			path = "../../../../Desktop/ToriBigData/rasters/"
# 		) ## Soil download not working, use VSI to interface in the meantime, read to local variable to resample
# 		# soil[[p]] <- soilworld
# 	
# 	
# }
# 
# ## Create a list for sets of variables downloaded with same naming convention, update as you go
# wc_var <- c("bio","srad","vapr")
# 
# ## Download for WorldClim Data
# 
# ## Create Empty list
# clim <- c()
# 
# ## Loop for WorldClim Data
# for (q in 1:length(wc_var)) {
# 	var <- wc_var[q]
# 	print(var)
# 	
# 	worldclim <- geodata::worldclim_global(var = var,
# 																				 res = .5,
# 																				 path = "../../../../Desktop/ToriBigData/rasters/")
# 	# clim[[q]] <- worldclim
# 	
# }



# ## Unlist and Stack Predictor rasters, use numerals to not overwrite OG rasters (Use Raster package if crop too large for terra)
# soil2 <- stack(rast(soil))
# clim2 <- stack(rast(unlist(clim)))
# 
# ## Extend Rasters to Match Largest Extent (soilgridsv2, in initial case)
# clim3 <- raster::crop(clim2,extent(soil2))
# 
# ## Stack rasters together with vector
# predictors <- c(rast(soil2),rast(clim3))


## AS OF APRIL 30 2024, THIS SCRIPT WON'T WORK USING ALL RASTERS, IF THIS IS TRUE FOR YOU, MOVE TO NEXT SCRIPT AND A SMALLER LIST OF RASTERS WILL BE GENERATED
rast_files <- list.files(paste0("../../../../Desktop/ToriBigData/rasters/"), pattern = "*.tif", full.names = TRUE)

## Stack Your Cropped Area Rasters, crop (by smaller extent) if extents differ
one <- raster::stack(rast_files[1:27])
two <- raster::crop(raster::stack(rast_files[28:71]), extent(one))

predictors <- stack(one, two)


## Next step, extract data points
## Crop to disinct regions,

s_ext <- s_datum %>% 
	st_as_sf(coords = c("x","y"), crs = 4326) %>% 
	st_bbox()

n_ext <- n_datum %>% 
	st_as_sf(coords = c("x","y"), crs = 4326) %>% 
	st_bbox()

## Crop to Extent of S and N Populations, Output = Raster Stack
s_ENVR <- crop(predictors, s_ext)
writeRaster(s_ENVR,"../../../../Desktop/ToriBigData/sENVR.grd", overwrite=T)
n_ENVR <- crop(predictors, n_ext)
writeRaster(n_ENVR,"../../../../Desktop/ToriBigData/nENVR.grd", overwrite=T)


## Exttract OOCs
s_extract <- cbind(s_datum, extract(s_ENVR, s_datum))
s_extract2 <- s_extract[,colSums(is.na(s_extract))<nrow(s_extract)]
s_extract2 <- s_extract2[complete.cases(data.frame(s_extract2)), ]

n_extract <- cbind(n_datum, extract(n_ENVR, n_datum))
n_extract2 <- n_extract[,colSums(is.na(n_extract))<nrow(n_extract)]
n_extract2 <- n_extract2[complete.cases(data.frame(n_extract2)), ]

## Generate a Matrix of Cropped Raster Values
s_ENVM <- getValues(stack(s_ENVR))
n_ENVM <- getValues(stack(n_ENVR))

# Complete Cases, Save
s_ENVM <- s_ENVM[complete.cases(s_ENVM), ]
saveRDS(s_ENVM, file= "../../../../Desktop/ToriBigData/sENVM.rds")
n_ENVM <- n_ENVM[complete.cases(n_ENVM), ]
saveRDS(n_ENVM, file= "../../../../Desktop/ToriBigData/nENVM.rds")

# produce global environmental background data
globalEnvM <- rbind(s_ENVM, n_ENVM)
saveRDS(globalEnvM, file= "../../../../Desktop/ToriBigData/globeENVM.rds")

  ## Generate Basic Statistics of both Populations
Sstats <- apply(s_datum, 2, function(x) c(min = min(x), median = median(x), mean = mean(x), max = max(x), sd = sd(x)))
Nstats <- apply(n_datum, 2, function(x) c(min = min(x), median = median(x), mean = mean(x), max = max(x), sd = sd(x)))

# combine both statistics for direct comparison, write to results folder
ENvStats <- rbind(Sstats, Nstats)
write.csv(ENvStats, "rst/ENVstats.csv", row.names = F)

## Spatial Autocorrelation
#S population
pdf("rst/Spop_correlogram.pdf")
ecospat.mantel.correlogram(dfvar=s_extract2[c(1:72)],colxy=1:2, n=100, colvar=3:72, 
													 max=10, nclass=10, nperm=100)
dev.off()
#N Population
pdf("rst/Npop_correlogram.pdf")
ecospat.mantel.correlogram(dfvar=n_extract2[c(1:72)],colxy=1:2, n=100, colvar=3:72, 
													 max=10, nclass=10, nperm=100)
dev.off()


# Modeling ----------------------------------------------------------------


pca.clim <- dudi.pca(globalEnvM, center = TRUE,
										 scale = TRUE, scannf = FALSE, nf = 2)

global.scores <- pca.clim$li

SLS.scores <-
	suprow(pca.clim,
				 data.frame(s_datum)[, colnames(globalEnvM)])$li   
NLS.scores <-
	suprow(pca.clim,
				 data.frame(n_datum)[, colnames(globalEnvM)])$li

SEnv.scores <- suprow(pca.clim, s_ENVM)$li
NEnv.scores <- suprow(pca.clim, n_ENVM)$li


data.frame(s_datum)[, colnames(globalEnvM)]
# calculate the Occurrence Density Grid for both native and invasive species
SpopGrid <- ecospat.grid.clim.dyn(global.scores,
																		SEnv.scores,
																		SLS.scores)

NpopGrid <- ecospat.grid.clim.dyn(global.scores,
																			NEnv.scores, 
																			NLS.scores)

## Plot Occurrence Density grid
ecospat.plot.niche.dyn(SpopGrid, NpopGrid, quant = 0.1, interest = 2, title = "Niche Overlap", name.axis1 = "PC1", name.axis2 = "PC2")



## List of Figures
## Histogram: Similarity vs Equivalency
## Niche Dynamics Centriods: Plot against most important variables from MAXENT

