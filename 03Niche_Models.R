# 03Niche_Models.R
# Step 3

# Tori Ford
# April 2, 2024
# Objective: Complete All Steps of Ecological NIche Models, From Maxent
# https://rpubs.com/janikh99/1044769
# Modified from tutorial, native dynamics vs invaded


# Load in relevant Libraries
library(ecospat)
library(geodata)
library(raster)
library(rangeBuilder)
library(maptools)
library(sf)
library(dplyr)
library(devtools)
library(terra)
library(ade4)
library(ape)
library(biomod2)



# Premodel ----------------------------------------------------------------
# During premodeling, you must organize both rasters and occurrences!

datum <- read.csv("raw/CJG_popk2_for_TF.csv")
datum <- datum[,c(3,4)] ## Subset out long~lat coordinates, whole dataset
s_datum <- read.csv("raw/Senv1.csv")
s_datum <- s_datum[,c(1,2)] ## Subset out long~lat coordinates, population dataset
n_datum <- read.csv("raw/Nenv2.csv")
n_datum <- n_datum[,c(1,2)] ## Subset out long~lat coordinates, population dataset


# Read in total area rasters
rast_files <- list.files(paste0("raw/rast/total"), pattern = "*.asc", full.names = TRUE)

# ## if extents differ, crop (by smaller extent) 
# one <- rast(rast_files[1:27])
# two <- terra::crop(rast(rast_files[28:47]), ext(one))

## Render raster in, use <- c() if rasters are already loaded
predictors <- rast(rast_files)

## Extract data from predictor layers per occ point
datum.occ <- cbind(datum, terra::extract(predictors, datum))
datum.occ <- datum.occ[complete.cases(data.frame(datum.occ)), ] # use complete cases to remove

## Extract per population
s.occ <- cbind(s_datum, terra::extract(predictors, s_datum))
s.occ <- s.occ[complete.cases(data.frame(s.occ)), ]
s.occ <- s.occ[,c(4:73)] # .occ's will only contain predictor data
n.occ <- cbind(n_datum, terra::extract(predictors, n_datum))
n.occ <- n.occ[complete.cases(data.frame(n.occ)), ]
n.occ <- n.occ[,c(4:73)] # .occ's will only contain predictor data


## Next step, calculate the extents of the populations

s_ext <- s_datum %>% 
	st_as_sf(coords = c("x","y"), crs = 4326) %>% 
	st_bbox()

n_ext <- n_datum %>% 
	st_as_sf(coords = c("x","y"), crs = 4326) %>% 
	st_bbox()


## Crop to Extent of S and N Populations, use clusters for faster processing

beginCluster(10)
s_ENVR <- crop(predictors, s_ext)
n_ENVR <- crop(predictors, n_ext)
endCluster()

## Extract values normally:
beginCluster(10)
s_ENVM <- values(s_ENVR)
n_ENVM <- values(n_ENVR)
endCluster()

## Write out rasters, this may take some time
writeRaster(s_ENVR,"rst/sENVR.grd", overwrite=T)
writeRaster(n_ENVR,"rst/nENVR.grd", overwrite=T)

# Complete Cases, Save
s_ENVM <- s_ENVM[complete.cases(s_ENVM), ]
saveRDS(s_ENVM, file= "raw/sENVM.rds")
n_ENVM <- n_ENVM[complete.cases(n_ENVM), ]
saveRDS(n_ENVM, file= "raw/nENVM.rds")

# produce global environmental background data
globalEnvM <- rbind(s_ENVM, n_ENVM)
saveRDS(globalEnvM, file= "raw/globeENVM.rds")

## Generate Basic Statistics of both Populations
Sstats <- apply(s_datum, 2, function(x) c(min = min(x), median = median(x), mean = mean(x), max = max(x), sd = sd(x)))
Nstats <- apply(n_datum, 2, function(x) c(min = min(x), median = median(x), mean = mean(x), max = max(x), sd = sd(x)))

# combine both statistics for direct comparison, write to results folder
ENvStats <- rbind(Sstats, Nstats)
write.csv(ENvStats, "rst/ENVstats.csv", row.names = F)

## Spatial Autocorrelation, if the points are heavily correlated, you may want to consider variable reducing methods. 
#S population
pdf("rst/Spop_correlogram.pdf")
ecospat.mantel.correlogram(dfvar=s_extract[c(1:48)],colxy=1:2, n=100, colvar=3:48, 
													 max=10, nclass=10, nperm=100)
dev.off()
#N Population
pdf("rst/Npop_correlogram.pdf")
ecospat.mantel.correlogram(dfvar=n_extract[c(1:48)],colxy=1:2, n=100, colvar=3:48, 
													 max=10, nclass=10, nperm=100)
dev.off()


# Modeling ----------------------------------------------------------------

# Run PCA on global value matrix
pca.clim <- dudi.pca(globalEnvM, center = TRUE,
										 scale = TRUE, scannf = FALSE, nf = 2)

global.scores <- pca.clim$li

# Find supplementary rows in PCA by matching to columns from occ dataset
SLS.scores <-
	suprow(pca.clim, s.occ[, colnames(globalEnvM)])$li   

NLS.scores <-
	suprow(pca.clim, n.occ[, colnames(globalEnvM)])$li

# Find supplementary rows in PCA by matching to columns from envm dataset
SEnv.scores <- suprow(pca.clim, s_ENVM)$li
NEnv.scores <- suprow(pca.clim, n_ENVM)$li



# calculate the Occurrence Density Grid for both native and invasive species
SpopGrid <- ecospat.grid.clim.dyn(global.scores,
																	SEnv.scores,
																	SLS.scores)

NpopGrid <- ecospat.grid.clim.dyn(global.scores,
																	NEnv.scores, 
																	NLS.scores)

## Plot Occurrence Density grid
pdf("rst/ESpace_Niche_Overlap.pdf")
ecospat.plot.niche.dyn(SpopGrid, NpopGrid, quant = 0.1, interest = 2, title = "Niche Overlap", name.axis1 = "PC1", name.axis2 = "PC2")
dev.off()

# plot variable contributions

ecospat.plot.contrib(contrib=pca.clim$co, eigen=pca.clim$eig)

# Expansion = sp2
# Intersection  = If intersection=NA, the analysis is performed on the whole environmental extent (native and invaded). If intersection=0, the analysis is performed at the intersection between native and invaded range. If intersection=0.05, the analysis is performed at the intersection of the 5th quantile of both native and invaded environmental densities.
ecospat.niche.dyn.index(SpopGrid, NpopGrid, intersection = NA)$dynamic.index.w


## Geographic projection of the selected "index" stats
geoProj <- ecospat.niche.dynIndexProjGeo(SpopGrid,
																				 NpopGrid,
																				 env = raster::stack(predictors),index="unfilling")


## G-space Projection
# Here, we define the bounding box of our area of interest. This will change depending on the populations of observation
geoGrid <- expand.grid(longitude =
											 	seq(20, 160, length.out = 250),
											 latitude =
											 	seq(-20, 80, length.out = 250))

# Open map dataset from "maps", subset a region out
data("wrld_simpl")
mask <- subset(wrld_simpl, REGION = c("142","150")) # Region for Central Asia

## Use the bounding boxes as grids, and dataset as cordinates, and countries for mask
s.GeoGrid <- ecospat.grid.clim.dyn(geoGrid, geoGrid,
																			coordinates(s_datum),
																			geomask = mask)
n.GeoGrid <- ecospat.grid.clim.dyn(geoGrid, geoGrid,
																	 coordinates(n_datum),
																	 geomask = mask)

# Plot and save as pdf file
pdf("rst/GSpace_Niche_Overlap.pdf")
ecospat.plot.niche.dyn(s.GeoGrid, n.GeoGrid, quant = 0)
plot(wrld_simpl, add = TRUE)
dev.off()


# Sim-Eq Tests ------------------------------------------------------------

# perform the Niche Equivalency Test

eq.test <- ecospat.niche.equivalency.test(SpopGrid, NpopGrid, rep = 100, ncores = 2)

# perform the Niche Similarity Test

sim.test <- ecospat.niche.similarity.test(SpopGrid, NpopGrid, rep = 100, rand.type = 2, ncores = 2)

# plot Equivalency and Similarity Test
pdf("rst/Niche_EQandSIM_NSPop.pdf")
par(mfrow=c(1,2))
ecospat.plot.overlap.test(eq.test, "D", "Equivalency") 
ecospat.plot.overlap.test(sim.test, "D", "Similarity")
dev.off()



# Dynamics ----------------------------------------------------------------

# Create a niche overlap directory for all response variable dynamics plots
dir.create("rst/NicheOverlapDyn")

## Maybe loop
for (x in 1:ncol(s.occ)) {
	print(colnames(s.occ)[x])
	
	# gridding the native niche
	grid.clim.t.s <- ecospat.grid.clim.dyn(glob = globalEnvM[,x],
																				 glob1 = data.frame(s_ENVM[,x]),
																				 (s.occ)[,x], R = 1000, th.sp = 0)
	
	# gridding the invasive niche
	grid.clim.t.n <- ecospat.grid.clim.dyn (glob = globalEnvM[,x], 
																					glob1 = data.frame(n_ENVM[,x]), 
																					(n.occ)[,x], R = 1000, th.sp = 0)
	
	
	pdf(paste0("rst/NicheOverlapDyn/",colnames(s.occ)[x],"_NiOv.pdf"))
	
	ecospat.plot.niche.dyn(grid.clim.t.s, grid.clim.t.n, quant=0.1, interest=2, title= "Niche Overlap", name.axis1=paste0("Average ", colnames(s.occ)[x]))
	
	# showing the shift of the niche centroid along the temperature gradient (compared to the shift of the available climate in the study area)
	ecospat.shift.centroids(data.frame(s.occ)[,x],
													data.frame(n.occ)[,x],
													data.frame(s_ENVM)[,x],
													data.frame(n_ENVM)[,x])
	dev.off()
	
}
