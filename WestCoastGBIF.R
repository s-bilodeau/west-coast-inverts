#clear environment
rm (list=ls())

#set working directory
setwd("")

#load packages
library(tidyverse)
library(rgbif)
library(terra)
library(sf)

#current data: GBIF.org (16 November 2026) GBIF Occurrence Download https://doi.org/10.15468/dl.m6dcbk
#this version uses Catologue of Life taxonomy
#details can be viewed at the DOI linked above
#the same dataset using the GBIF Backbone taxonomy (16 November 2026, https://doi.org/10.15468/dl.928jhp) 
#can be accessed by swapping the id number below to '0011061-261110092616502'

#NOTE: we are pulling occurrence data from a larger range than our intended study area 
#in order to distinguish between actual range edges and study area edges later on

#download and unzip the gbif dataset
occurrences <- occ_download_get('0010962-261110092616502') %>%
  occ_download_import()
#occurrences <- occ_download_import(key="0010962-261110092616502", path=".") #use this line if offline and the zip file is in the working directory already

#check total species
length(unique(occurrences$taxonKey))


###for each observation, get depth at that point from the bathymetry###

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

#clean up the environment by removing unnecessary objects
#rm(depths, depth) #only necessary if memory is limited

#there are several rows with NA depth values...let's find out why!
nodepth <- occ[which(is.na(occ$ETOPO2022NEPacificIceSurfaceArcMinute)==TRUE),]
#all appear to be on the Colombia/Panama border and can be safely omitted 
#this is likely due to a mismatch in cropping our raster and our GBIF extent  
#but doesn't matter, since both go well outside our intended final extent
occ <- occ[which(is.na(occ$ETOPO2022NEPacificIceSurfaceArcMinute)==FALSE),]

#there are a lot more rows missing GBIF depth
nodepth <- occ[which(is.na(occ$depth)==TRUE&is.na(occ$elevation)==TRUE),]
#in theory, any value over zero is on land, but there are likely intertidal species with positive depth/elevation values due to poor location mapping

#we want to eliminate observations too far north/south to fall in our study area (with a little wiggle room), 
#ideally before grouping by spp so that spp with no observations in the study area are excluded
occ_short <- occ[which(occ$decimalLatitude > 22 & occ$decimalLatitude < 48.4),] #this cuts down the total unique taxonKey (spp) values

#clean up the environment by removing unnecessary objects
#rm(occ, nodepth) #only necessary if memory is limited

#multiply ETOPO (raster) depth column and GBIF elevation column by -1, so depth values are positive and elevation values negative (to match with GBIF reported depths)
occ_short$ETOPO2022NEPacificIceSurfaceArcMinute <- (-1)*as.numeric(occ_short$ETOPO2022NEPacificIceSurfaceArcMinute)
occ_short$elevation <- (-1)*as.numeric(occ_short$elevation)

#we will convert infinite depth values to NAs
occ_short$depth[is.infinite(occ_short$depth)] <- NA

#fill missing (NA) GBIF depth data from the GBIF elevation column (first choice) or the (transformed) ETOPO depth column; this takes a while to run!
for(i in 1:nrow(occ_short)) {
  if(is.na(occ_short[i,"depth"])) { 
    if(is.na(occ_short[i,"elevation"])){
    occ_short[i,"depth"] <- as.numeric(occ_short[i,"ETOPO2022NEPacificIceSurfaceArcMinute"]) 
    }
    else { 
    occ_short[i,"depth"] <- as.numeric(occ_short[i,"elevation"]) 
    }
  }
}

#save the output after this for loop since it takes a *long* time to run!
saveRDS(occ_short, "cropped_occ_with_depth.rds")

#occ_short <- readRDS("cropped_occ_with_depth.rds") #un-comment this line to read in the saved data, 
#rather than (re)running the time-consuming code up to this point
#the .RDS file is not included in the GitHub repo due to size, so you will have to generate it yourself the first time


###now, we want to eliminate species with fewer than 25 observations and flag those with extreme depth/elevation values###

#first, create a column to count total observations for each species (taxonKey), 
#then subset the dataframe to only species with 25+ observations
occ_short$obs <- as.numeric(ave(as.character(occ_short$taxonKey), as.character(occ_short$taxonKey), FUN=length))
occ_25 <- occ_short[which(occ_short$obs>=25),]

#we also want to drop taxa that are not identified to at least the genus level (i.e., genus is empty "")
occ_25 <- occ_25[-which(occ_25$genus==""),]

#now, we want to group observations by species (taxonKey) 
#and count how many of those observations have depth >200m or <0 (i.e., above sea level)
deep <- occ_25 %>%
  group_by(taxonKey, scientificName, phylum, class, obs) %>%
  summarize(count=sum(depth>200))
deep <- deep[which(deep$count>(as.numeric(deep$obs)/4)),] #keep a record of only species for which >25% (1/4) obs are >200m deep
  
land <- occ_25 %>%
  group_by(taxonKey, scientificName, phylum, class, obs) %>%
  summarize(count=sum(depth<0))
land <- land[which(land$count>(as.numeric(land$obs)/2)),] #keep a record of only species for which >50% (1/2) obs have negative depths 
#(i.e., are above sea level)

###we will now manually check the flagged land group and eliminate all species from the deep group entirely###

#write.csv(land, "potentialTerrestrialSpeciesToCheck.csv", row.names=TRUE) #use to write a copy of the land data for manual checking

#after manually checking all species, read in a copy of the land data (with an added binary indicator column "marine" (1=yes, 0=no))
land <- read.csv("potentialTerrestrialSpeciesToCheck.csv") #data is sorted by phylum, class, and scientificName; X is the original row number before sorting

#notes from manual filtering:
#brackish/fresh species are flagged as non-marine unless additional info supports their retention
#brackish only and brackish/marine species are retained as marine in the dataset
#genera with both marine and freshwater where the species is unknown were flagged as marine for retention in the absence of other details
#Uncinocythere spp. are ostracods that live on crayfish and are presumed to be non-marine
#Aratus pisonii is the mangrove tree crab (likely relevant since larvae are marine and planktonic after brooding of eggs, either misidentified or introduced on the West Coast)
#Cardisoma crassum is another terrestrial/marine crab
#Coenobita compressus is a terrestrial/marine hermit crab
#Pohlia appears to be a fossil genus despite fossils being filtered out on GBIF and was therefore flagged for exclusion
#many nematodes were only identified to the genus level (more than in other taxonomic groups) and many of these genera are both freshwater and marine

land <- land[which(land$marine==1),] #drop non-marine species

#filter overall dataset based on depth and terrestrial filtering steps
occ_25_filtered <-  occ_25[occ_25$taxonKey %in% land$taxonKey,] #keeping only marine species identified above
occ_25_filtered <-  occ_25[!occ_25$taxonKey %in% deep$taxonKey,] #dropping all species with >25% occurrences below 200m

#write an updated .csv of the remaining species after all our filtering
spp25 <- occ_25_filtered[,c("phylum", "class", "order", "family", "genus", "species", "scientificName", "taxonKey")] %>% 
      group_by(taxonKey, scientificName, phylum, class, order, family, genus, species) %>%
      summarize(
        obs=n()
      )
#write.csv(spp25, "AllWestCoastSppSummarized.csv", row.names=FALSE)








############# Code Below is for getting depth and range values for all remaining species and is currently incomplete ############
#
####summarize observations by species, recording the min, max, and median depth for each species (by both measurement systems) and the number of observations###
#  
##(we will eventually also want north/south boundaries, but not until we further constrain our area of interest)
##NOTE: there are many rows where individualCount is a value <1, so we won't attempt to quantify observed individuals this way 
#  
##summarize with max depths
#spp <- occ_short[,c(5,6,13,28,34,51)] %>% #phylum, class, scientificName, individualCount, depth, taxonKey, ETOPO2022NEPacificIceSurfaceArcMinute (AKA depth from bathymetry)
#  group_by(taxonKey, scientificName, phylum, class) %>%
#  summarize(
#    across(matches("depth|ETOPO2022NEPacificIceSurfaceArcMinute"), max, na.rm=T),
#    obs=n()
#  )
##summarize with min depths
#sppmin <- occ_short[,c(5,6,13,28,34,51)] %>% #phylum, class, scientificName, depth, taxonKey, ETOPO2022NEPacificIceSurfaceArcMinute (AKA depth from bathymetry)
#  group_by(taxonKey, scientificName, phylum, class) %>%
#  summarize(
#    across(matches("depth|ETOPO2022NEPacificIceSurfaceArcMinute"), min, na.rm=T),
#    obs=n()
#  )
##summarize with median depths
#sppmed <- occ_short[,c(5,6,13,28,34,51)] %>% #phylum, class, scientificName, depth, taxonKey, ETOPO2022NEPacificIceSurfaceArcMinute (AKA depth from bathymetry)
#  group_by(taxonKey, scientificName, phylum, class) %>%
#  summarize(
#    across(matches("depth|ETOPO2022NEPacificIceSurfaceArcMinute"), median, na.rm=T),
#    obs=n()
#  )
##rename columns and join max and min depth summaries
#names(spp) <- c("taxonKey",  "scientificName", "phylum", "class", "GBIF_depth_max", "ETOPO_depth_max", "Observations")
#names(sppmin) <- c("taxonKey",  "scientificName", "phylum", "class", "GBIF_depth_min", "ETOPO_depth_min", "Observations")
#names(sppmed) <- c("taxonKey",  "scientificName", "phylum", "class", "GBIF_depth_med", "ETOPO_depth_med", "Observations")
#spp <- left_join(spp, sppmin)
#spp <- left_join(spp, sppmed)
#
#
#
##our summarize function may have re-introduced some Inf values, despite filtering for them earlier, so let's remove those!
#spp25$GBIF_depth_min[is.infinite(spp25$GBIF_depth_min)] <- NA 
#spp25$GBIF_depth_med[is.infinite(spp25$GBIF_depth_med)] <- NA    
#spp25$GBIF_depth_max[is.infinite(spp25$GBIF_depth_max)] <- NA    
#
#
#hist(spp25$GBIF_depth_med) #most remaining species have a median depth of 500m or less
#shallow <- spp25[which(spp25$GBIF_depth_med<50),] 
#hist(shallow$GBIF_depth_med) #while many shallower spp have median depth <20m, 
##there is no clear bimodal distribution like Pappalardo et al. describe in their data





