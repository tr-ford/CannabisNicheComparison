# ENMEval.R
# Step 2

# Tori Ford
# April 3, 2024

library(ENMeval)
library(raster)
library(biomod2)
library(dismo) 
library(ecospat)
## Optional Java Options, if you are encountering memory issues, this must be changed per your needs: RUN B4 rJava, if not do rm(list = ls()), then try again.
options(java.parameters = "-Xmx64g")
library(rJava) 
library(rasterVis)
library(latticeExtra)
library(parallel) 
library(sf)
library(terra)
library(dplyr)
library(sp)
library(gtools)
library(usdm)

## Dataframe of ALL K groups, rename columns per need
ssp <- read.csv("raw/CJG_popK2_for_TF.csv")
ssp <- ssp %>% rename(x = LON,y = LAT)
ssp$species <- "Cannabis sativa"

## Read in Rasters of Shared Region
# rast_files <- list.files(paste0("raw/rast/buf/"), pattern = "*.asc", full.names = TRUE)
rast_files <- list.files(paste0("../../../../Desktop/ToriBigData/rasters/"), pattern = "*.tif", full.names = TRUE)
totalr_files <- list.files(paste0("raw/rast/total/"), pattern = "*.asc", full.names = TRUE)

## Stack Your Cropped Area Rasters, crop (by smaller extent) if extents differ
soils <- raster::stack(rast_files[1:27])
clims <- raster::crop(raster::stack(rast_files[28:71]), extent(soils))

envs <- stack(soils, clims)

## Define Your Projection String
projection(envs) <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs" 

## IF YOU NEED TO CLEAR RASTER TEMPS, USE removeTmpFiles(h=0)
set.seed(333) # Testing reproducibility with set seed

for (s in unique(ssp$POPk2)) {
	
	## Set up buffer for crop + mask
	## Total Extent
	ssp.sf <- sf::st_as_sf(ssp, coords = c("x","y"), crs = raster::crs(envs))
	
	
	## Get bbox for raster crop to K group
	ssp.ext <- ssp.sf %>% 
		st_as_sf(coords = c("x","y"), crs = 4326) %>% 
		st_bbox()
	
	## Create Raster of Total Extent, For Projection Step
	## Check for existing total extent raster files, no need to crop if they exist
	if (length(rast_files) == length(totalr_files)) {
		bbox_envs <- stack(totalr_files)
		print(Sys.time())
	} else {
 		bbox_envs <- raster::crop(envs, extent(ssp.ext))
		bbox_envs <- raster::crop(bbox_envs, extent(ssp.ext))
		for (b in 1:length(names(bbox_envs))) {
			layer <- bbox_envs[[b]]
			name <- names(layer)
			raster::writeRaster(layer, paste0("../../../../Desktop/ToriBigData/rasters/total/",name,".asc"), overwrite = T)
			# raster::writeRaster(layer, paste0("raw/rast/total/",name,".asc"), overwrite = T)
			print(Sys.time())
		}
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


	buffers <- c()
	dir.create(paste0("../../../../Desktop/ToriBigData/rasters/buf",s))
	# dir.create(paste0("raw/rast/buf",s))
	## CROP, MASK, SAVE RADIAL BUFFER
	for (o in 1:length(rast_files)) {
		# Subset raster layer
		layer <- rast(envs[[o]])
		# Setup file names
		name <- names(layer)
		print(name) ## Just here to check on progress of loop
		# out <- paste0(path, name)
		# outfile <- paste0(out, end)
		# Crop and mask
		cropped <- terra::crop(layer, xy.buf) #ext
		masked <- terra::mask(cropped, xy.buf)
		# buffers[[o]] <- raster(masked)
		# Write raster
		# terra::writeRaster(masked,
		# 									 filename = paste0("raw/rast/buf",s,"/",name,".asc"),
		# 									 overwrite = F)
		terra::writeRaster(masked,
											 filename = paste0("../../../../Desktop/ToriBigData/rasters/buf",s,"/",name,".asc"),
											 overwrite = F)
		print(paste0(o))
	}
	
	
		 ## IF NEEDED, ADD VIF (OR OTHER VARIABLE REDUCING STEPS) HERE
	
	## Stack Radial Buffers
	# radial.buffers <- stack(list.files(paste0("raw/rast/buf",s,"/"), pattern = "*.asc", full.names = TRUE))
	radial.buffers <- stack(list.files(paste0("../../../../Desktop/ToriBigData/rasters/buf",s,"/"), pattern = "*.asc", full.names = TRUE))
	## CHECK ENMEVAL VERSION
	
	## All of the following steps were derived from the USDM Vignette: https://cran.r-project.org/web/packages/usdm/usdm.pdf, see 'vif' and 'exclude'.
	#calculate VIFs using a threshold of 10
	vif.buffers <- vifstep(radial.buffers, th=10) #USDM Vignette: "A VIF greater than 10 is a signal that the model has a collinearity problem".
	print("step complete")
	
	## Get result table with VIF scores for kept variables
	vif.table <- as.data.frame(vif.buffers@results)
	write.csv(vif.table, paste0("rst/",s,"pop_VIF_Results.csv"))
	
	## Use the 'exclude' function to exclude rasters past threshold and stack the remainder. 
	in.vif <- exclude(radial.buffers,vif.buffers)
	print("layers excluded")
	
	## Also exclude nodata values
	NaN_DF <- vif.buffers@results[!is.na(vif.buffers@results$VIF),]
	maxent.vif <- raster::subset(in.vif,NaN_DF$Variables)
	
	## Randomize background points
	bg <- dismo::randomPoints(maxent.vif[[1]], n = 10000) %>% as.data.frame() 
	colnames(bg) <- colnames(xy)
	
	## Create jacknife partition by sampling background points against occurrences in lon_lat format
	checkr <- get.checkerboard1(xy,maxent.vif, bg, aggregation.factor = 4)
	table(checkr$occs.grp)
	
	EcologicalNiche <- ENMevaluate(occs = xy, 
													envs = maxent.vif, 
													bg = bg, 
													algorithm = "maxent.jar", 
													partitions = "checkerboard1", 
													tune.args = list(fc = c("L","LQ","LQH","H","QH","Q","LH"),
																					 rm = 1:5),
													parallel = TRUE,
													numCores = 36)
	## Save your Maxent Model
	saveRDS(EcologicalNiche, file = paste0("raw/", s, "_ENMeval.RDS"))
	
	
	## Calculate Niche Overlap between models,which support Schoener’s D, to compare output similarity
	overlap_D <- calc.niche.overlap(EcologicalNiche@predictions, overlapStat = "D")
	write.table(overlap_D, file = paste0("rst/", s, "_Schoeners_D_Niche_Overlap.txt"), sep = "\t")
	## Calculate Niche Overlap between models,which support Moran’s I, to compare output similarity
	overlap_I <- calc.niche.overlap(EcologicalNiche@predictions, overlapStat = "I")
	write.table(overlap_I, file = paste0("rst/", s, "_Morans_I_Niche_Overlap.txt"), sep = "\t")
	
	
	# Overall results Table
	results <- eval.results(EcologicalNiche)
	## Save model results
	write.table(results, file = paste0("rst/", s, "_Model_Results.txt"), sep = "\t")
	
	
	
	## Use AICc to filter for an optimal model
	aicc.rst <- results %>% dplyr::filter(delta.AICc == min(na.exclude(delta.AICc))) #%>% filter(delta.AICc == 0)
	write.table(aicc.rst, file = paste0("rst/", s, "_AICc_Model_Results.txt"), sep = "\t")
	
	## Create Raster of Selected Model, save it
	aicc.raster <- eval.predictions(EcologicalNiche)[[aicc.rst$tune.args]]
	writeRaster()
	
	## Select a single model using the parameters of our aicc.rst filter
	aicc.model <- eval.models(EcologicalNiche)[[aicc.rst$tune.args]]
	saveRDS(aicc.model, paste0(s,"_ENM.RDS"))
	plot(aicc.model, type = "cloglog")
	
	## Generate Null Models for Validation, Use get() to pull from results table
	null.models <- ENMnulls(e = EcologicalNiche, mod.settings = list(fc = get(aicc.rst$a), rm = get(aicc.rst$b)), no.iter = 200)
	
	## Extract the null.models results into a table
	null.table <- null.emp.results(null.models)
	write.table(null.table, file = paste0("rst/", s, "_Null_Results.txt"), sep = "\t")
	
	## Create figures, save as vectored graphics.
	svg(filename = paste0("rst/", s, "_Null_Figures.svg"), width = 12, height = 12)
	evalplot.nulls(null.models, stats = c("or.10p", "auc.val","auc.train","auc.diff","cbi.val"), plot.type = "histogram")
	dev.off()	
	
	
	
	## Convert the raster to spatial points
	raster.points <- rasterToPoints(aicc.raster, spatial = TRUE)
	
	## Convert spatial points to a dataframe
	raster.df  <- data.frame(raster.points)
	
	## Get Extent of raster to constrain map to region
	exts <- extent(radial.buffers)
	
	## Call world map (can restrict by region, refer to plotting scheme in 06_ENM_Evaluation.R)
	regmap <- map_data("world",xlim = c(exts@xmin,0),ylim = c(exts@ymin,exts@ymax))
	
	## Create a map of the optimal model, plotted on top of a world map and add points too represent occurrences
	model_map <- ggplot() +
		geom_polygon(data = regmap, aes(x = long, y = lat, group = group),color = "snow", fill = "grey38")+
		geom_raster(data = raster.df , aes(x = x, y = y, fill = layer, alpha = 1)) + 
		scale_fill_viridis(option = "G",direction = 1, name = "Habitat Suitability") +
		scale_alpha(guide = 'none')+
		geom_point(data=ssp, aes(x=x, y=y, color = species),shape=8,stroke = 0.4,size = 2,alpha = 1, show.legend = TRUE) + 
		scale_color_manual(values = c('#355A20')) +
		ggtitle(paste0(s, " Population Subset of Cannabis sativa Habitat Suitability"))+ 
		coord_quickmap()
	
	## Save the globally/regionally projected model (change name per needs)
	ggsave(paste0("ENM_Figure_for_",s,"_Subset.svg"),	plot = model_map ,	path = "rst/",	height = 8,	width = 12)
	
	
	
	## Use DISMO to project the model into the total area between populations
	total_proj <- dismo::predict(aicc.model, bbox_envs, progress = 'text', overwrite = TRUE,filename = paste0("rst/",s, "_Total_Projection.grd"))#
	
	print(paste0(Sys.time()," Loop Done  "))
	}
  