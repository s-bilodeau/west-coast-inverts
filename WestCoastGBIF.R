#clear environment
rm (list=ls())

#set working directory
setwd("~/WestCoastSpecies/west-coast-inverts")

#load packages
library(tidyverse)
library(rgbif)
library(terra)
library(sf)

#current data: GBIF.org (11 November 2026) GBIF Occurrence Download https://doi.org/10.15468/dl.hd6k4u
#this version uses Catologue of Life

#unzip and read in gbif download
occurences <- occ_download_get('0002816-261110092616502') %>%
  occ_download_import()
#occurrences <- occ_download_import(key="0002816-261110092616502", path=".") #use this line if the zip file is in the working directory already

#check total species
length(unique(occurrences$taxonKey))

#for each observation, get depth at that point from the bathymetry

#make a spatial object with just decimalLatitude and decimal Longitude and gbifID (unique to each row/observation)
pts <- st_as_sf(occurrences[,c(1,22,23)], coords=c("decimalLongitude","decimalLatitude"), crs=4326)

#read in raster with depth/elevation for the northeast Pacific (Panama up to Alaska)
#this is from the ETOPO 2022 Arc Minute Scale Ice Surface Topology model
r <- rast("ETOPO2022NEPacificIceSurfaceArcMinute.tiff", drivers="GTiff") #GeoTIFF from the Ice Surface version of the 2022 ETOPO 60 arc second (1 arc minute) model
#this raster's CRS is latlon WGS84, EPSG 4326
#GBIF data/coordinates uses this same CRS (EPSG 4326) by default, which is very convenient!

#extract depth values for each point location on the raster
depths <- extract(r, pts) #NOTE: if this fails, check that tidyverse packages that mask extract() weren't loaded after the terra library

#re-combine depth values with coordinates and gbifIDs
depths <- cbind(depths, pts)

#take just depth and gbifID and tack it onto the original occurrences dataframe
depth <- as.data.frame(depths[,2:3])
occ <- full_join(occurrences, depth, by="gbifID")

#there are apparently 13 rows with NA depth values...let's find out why!
nodepth <- occ[which(is.na(occ$ETOPO2022NEPacificIceSurfaceArcMinute)==TRUE),]
#all appear to be on the Colombia/Panama border and can be safely omitted 
#this is likely due to a mismatch in cropping our raster and our GBIF extent but doesn't really matter, 
#since both go well outside our intended final extent
occ <- occ[which(is.na(occ$ETOPO2022NEPacificIceSurfaceArcMinute)==FALSE),]

#there are a lot more rows missing GBIF depth and/or elevation
nodepth <- occ[which(is.na(occ$depth)==TRUE&is.na(occ$elevation)==TRUE),]
#in theory, any value over zero is likely on land, but there are likely intertidal species with positive depth/elevation values due to poor location mapping

#we want to eliminate observations too far north/south to fall in our study area (with a little wiggle room), 
#ideally before grouping by spp so that spp with no observations in the study area are excluded
occ_short <- occ[which(occ$decimalLatitude > 22 & occ$decimalLatitude < 48.4),] #this cuts the total unique taxonKey (spp) values to 24023

#summarize observations by species, recording the min, max, and median depth for each species (by both measurement systems) and the number of observations
#(we will eventually also want north/south boundaries, but not until we further constrain our area of interest)
#NOTE: there are multiple (almost 2000) rows where individualCount is a value <1, so we won't attempt to quantify observed individuals this way 

#multiply ETOPO depth column by -1, so depth values are positive and elevation values negative (to match with GBIF reported depths)
occ_short$ETOPO2022NEPacificIceSurfaceArcMinute <- (-1)*occ_short$ETOPO2022NEPacificIceSurfaceArcMinute

#we will convert infinite depth values to NAs
occ_short$depth[is.infinite(occ_short$depth)] <- NA    

#summarize with max depths
spp <- occ_short[,c(5,6,13,28,34,51)] %>% #phylum, class, scientificName, individualCount, depth, taxonKey, ETOPO2022NEPacificIceSurfaceArcMinute (AKA depth from bathymetry)
  group_by(taxonKey, scientificName, phylum, class) %>%
  summarize(
    across(matches("depth|ETOPO2022NEPacificIceSurfaceArcMinute"), max, na.rm=T),
    obs=n()
  )
#summarize with min depths
sppmin <- occ_short[,c(5,6,13,28,34,51)] %>% #phylum, class, scientificName, depth, taxonKey, ETOPO2022NEPacificIceSurfaceArcMinute (AKA depth from bathymetry)
  group_by(taxonKey, scientificName, phylum, class) %>%
  summarize(
    across(matches("depth|ETOPO2022NEPacificIceSurfaceArcMinute"), min, na.rm=T),
    obs=n()
  )
#summarize with median depths
sppmed <- occ_short[,c(5,6,13,28,34,51)] %>% #phylum, class, scientificName, depth, taxonKey, ETOPO2022NEPacificIceSurfaceArcMinute (AKA depth from bathymetry)
  group_by(taxonKey, scientificName, phylum, class) %>%
  summarize(
    across(matches("depth|ETOPO2022NEPacificIceSurfaceArcMinute"), median, na.rm=T),
    obs=n()
  )
#rename columns and join max and min depth summaries
names(spp) <- c("taxonKey",  "scientificName", "phylum", "class", "GBIF_depth_max", "ETOPO_depth_max", "Observations")
names(sppmin) <- c("taxonKey",  "scientificName", "phylum", "class", "GBIF_depth_min", "ETOPO_depth_min", "Observations")
names(sppmed) <- c("taxonKey",  "scientificName", "phylum", "class", "GBIF_depth_med", "ETOPO_depth_med", "Observations")
spp <- left_join(spp, sppmin)
spp <- left_join(spp, sppmed)

#remove spp with fewer than 25 observations
spp25 <- spp[which(spp$Observations>24),] #this drops our total spp to a much more reasonable number!

#important note about depth: GBIF depths are in m below sea level, and ETOPO has been adjusted to match,
#so the max value is the deepest and the min the shallowest (including negative numbers indicating elevation above sea level)

#we want to fill in a "canonical" depth for those lacking GBIF depth/elevation data
#since GBIF is likely more accurate, we will use ETOPO for this only for entries otherwise lacking in depth/elevation data

#our summarize function may have re-introduced some Inf values, despite filtering for them earlier, so let's remove those!
spp25$GBIF_depth_min[is.infinite(spp25$GBIF_depth_min)] <- NA 
spp25$GBIF_depth_med[is.infinite(spp25$GBIF_depth_med)] <- NA    
spp25$GBIF_depth_max[is.infinite(spp25$GBIF_depth_max)] <- NA    

#copy GBIF depths to the canonical depth columns, then fill NA values with ETOPO data (transformed to be non-negative for depths)
spp25$depth_min <- spp25$GBIF_depth_min
spp25$depth_med <- spp25$GBIF_depth_med
spp25$depth_max <- spp25$GBIF_depth_max

for(i in 1:nrow(spp25)) {
  if(is.na(spp25[i,12])) {
    spp25[i,12] <- as.numeric(spp25[i,9])
  }
}

for(i in 1:nrow(spp25)) {
  if(is.na(spp25[i,13])){
    spp25[i,13] <- spp25[i,11]
  }
}

for(i in 1:nrow(spp25)) {
  if(is.na(spp25[i,14])){
    spp25[i,14] <- spp25[i,6]
  }
}


#we want to exclude spp with >25% of occurrences deeper than 200m
#we'll start by dropping any spp with a min depth (i.e., shallowest depth) >200m
spp25 <- spp25[which(spp25$depth_min<200),]

hist(spp25$GBIF_depth_med) #most remaining species have a median depth of 200m or less, so let's spot-check the deeper ones
shallow <- spp25[which(spp25$GBIF_depth_med<50),] 
hist(shallow$GBIF_depth_med) #while many spp have median depth <20m, there is no clear bimodal distribution like Pappalardo et al. describe in their data

#let's also flag any species with a negative maximum depth value, i.e., all observations are above sea level
spp25[which(spp25$depth_min<0),] #we don't appear to have any, which is good news!

#save as a .csv
#write.csv(spp25, "AllWestCoastSppSummarized.csv", row.names=FALSE)




