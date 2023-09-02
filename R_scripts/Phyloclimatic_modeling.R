#Phyloclimatic modeling
#Load the required packages for the analysis
library(raster)
library(phytools)
library(machuruku)
# Acquire environmental variables
bio01 = raster("D:/climate_data/bio01.tif")
bio07 = raster("D:/climate_data/bio07.tif")
bio12 = raster("D:/climate_data/bio12.tif")
bio15 = raster("D:/climate_data/bio15.tif")
Clim_Cur <- stack(bio01, bio07, bio12, bio15)
#Prepare presence locations
cza_occ <- read.delim("cza_occ.csv", h=T, sep=",")
#Load treefile
cza_tree <- read.nexus("cza_tree.treefile")
#Estimate tip response curves and ancestral niches (using 3.79 Mya as an example)
resp_cza <- machu.1.tip.resp(cza_occ , Clim_Cur, verbose = T)
ace_cza_t379 <- machu.2.ace(resp_cza, cza_tree, timeslice=3.79, unc=T)
#Projecte ancestral models into paleoclimatic data
# Acquire paleoclimatic data
#379 Mya
clim.t379.file <- list.files(path = "D:/Oscillayers/t379", full.names = TRUE)
Clim_t379 <- stack(clim.t379.file)
#Make projection
mod.t379 <- machu.3.anc.niche(ace_cza_t379, Clim_t379, verbose = T)
