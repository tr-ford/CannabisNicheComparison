# Trial_Invasion_Comparison.R
# Step 1

# Tori Ford
# Started March 19, 2024
# Modified from "humboldt" vignette: https://jasonleebrown.github.io/humboldt/

library(dplyr)
library(humboldt)
library(raster)
library(dismo)

datum <- read.csv("raw/CJG_popK2_for_TF.csv")

## Rename columns
datum <- datum %>% rename(
	x = LON,
	y = LAT
)

## Subset data into predetermined groups
s_datum <- datum[grep("S", datum$POPk2),]
n_datum <- datum[grep("N", datum$POPk2),]


## Combine Columns for Name Clarity (Gsub for Spaces)
s_datum$sp <- gsub(" ","", paste(paste0(s_datum$ID), "_", paste0(s_datum$POPk2)))
n_datum$sp <- gsub(" ","", paste(paste0(n_datum$ID), "_", paste0(n_datum$POPk2)))

## Generate Species Frames, Reaarange
s_sp1 <- s_datum[,c(3:4,79)]
s_sp1 <- dplyr::select(s_sp1,3,1,2)
write.csv(s_sp1, "raw/Ssp1.csv", row.names = F)
n_sp2 <- n_datum[,c(3:4,79)]
n_sp2 <- dplyr::select(n_sp2,3,1,2)
write.csv(n_sp2, "raw/Nsp2.csv", row.names = F)


## Generate ENVs Frames
s_env1 <- s_datum[,c(3:5,10:78)]
write.csv(s_env1, "raw/Senv1.csv", row.names = F)
n_env2 <- n_datum[,c(3:5,10:78)]
write.csv(n_env2, "raw/Senv2.csv", row.names = F)


## Scrub frames for NAs and false numbers
s_env1 <- humboldt.scrub.env(s_env1)
n_env2 <- humboldt.scrub.env(n_env2)


## Scrub species frames for NAs
s_sp1 <- na.exclude(s_sp1)
n_sp2 <- na.exclude(n_sp2)


## Convert to Numeric
s_env1$x <- as.numeric(s_env1$x)
s_sp1$x <- as.numeric(s_sp1$x)
s_env1$y <- as.numeric(s_env1$y)
s_sp1$y <- as.numeric(s_sp1$y)
s_sp1$sp <- as.factor(s_sp1$sp)

n_env2$x <- as.numeric(n_env2$x)
n_sp2$x <- as.numeric(n_sp2$x)
n_env2$y <- as.numeric(n_env2$y)
n_sp2$y <- as.numeric(n_sp2$y)
n_sp2$sp <- as.factor(n_sp2$sp)

## Reduce Variables for Top Contrib Model (Estimates?)
reduc.vars<- humboldt.top.env(env1 = s_env1, env2 = n_env2,
															sp1 = s_sp1, sp2 = n_sp2,
															rarefy.dist = 0,  
															env.reso = 0.00833333, learning.rt1 = 0.0001,
															learning.rt2 = 0.0001, e.var = (3:72),
															pa.ratio = 4, steps1 = 1, steps2 = 1,
															method="contrib", contrib.greater=4)


zz<-humboldt.g2e(env1 = s_env1, env2 = n_env2, 
								 sp1 = s_sp1, sp2 = n_sp2, 
								 reduce.env = 0, reductype = "PCA", 
								 non.analogous.environments = "YES", 
								 env.trim= T,env.trim.type = "RADIUS", e.var=c(3:72),  
								 col.env = e.var, trim.buffer.sp1 = 200, 
								 trim.buffer.sp2 = 200, rarefy.dist = 2,
								 rarefy.units = "km", env.reso = 0.00833333, 
								 kern.smooth = 1, R = 100, run.silent = F)


scores.env1<-zz$scores.env1[1:2]
scores.env2<-zz$scores.env2[1:2]
scores.env12<- rbind(zz$scores.env1[1:2],zz$scores.env2[1:2])
scores.sp1<-zz$scores.sp1[1:2]
scores.sp2<-zz$scores.sp2[1:2]

##estimate the Potential Niche Truncation Index
pnt1<- humboldt.pnt.index(scores.env12,scores.env1,scores.sp1,kern.smooth=1,R=100)
pnt2<- humboldt.pnt.index(scores.env12,scores.env2,scores.sp2,kern.smooth=1,R=100)

## run create a grid of Environmental Space Function
z1<- humboldt.grid.espace(scores.env12,scores.env1,scores.sp1,kern.smooth=1,R=100)
z2<- humboldt.grid.espace(scores.env12,scores.env2,scores.sp2,kern.smooth=1,R=100)

## plot niche in espace
a<- humboldt.plot.niche(z1,"Species 1","PC1","PC2")
humboldt.plot.niche(z2,"Species 2","PC1","PC2")

##run it first with full environmental for background tests and equivalence statistic (total equivalence or divergence in current distributions)
full<-humboldt.doitall(inname="full_extent", 
											 env1 = s_env1, env2 = n_env2, 
											 sp1 = s_sp1, sp2 = n_sp2, 
											 rarefy.dist=1, rarefy.units="km", 
											 env.reso=0.083333, reduce.env=0, 
											 reductype="PCA", non.analogous.environments="YES", 
											 correct.env=T, env.trim=F,  
											 pcx=1, pcy=2, 
											 col.env=e.var, e.var=c(3:72), 
											 R=100, kern.smooth=1, 
											 e.reps=200, b.reps=200, 
											 nae="YES",thresh.espace.z=0.0001, 
											 p.overlap=T, p.boxplot=T, 
											 p.scatter=T, run.silent=F, ncores=4)
