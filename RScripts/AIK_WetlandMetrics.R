#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Title: Wetland Metrics 
#Coder: C. Nathan Jones (cnjones7@ua.edu)
#Date: 6/12/2019
#Purpose: Develop a suite of spatial metrics for AIK's masters thesis
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#A few quick notes:
#1) All data is housed in the Palmer Lab subdirectory.  (Choptank/Nate/Storage_Capacity)
#2) This analysis relies on WBT GIS platform. [https://www.uoguelph.ca/~hydrogeo/Whitebox/download.shtml
#     Before running this script, users will need to download WBT and define the 
#     location of both the WBT executable and a scratch workspace

#1.0 Setup workspace~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Clear memory
rm(list=ls(all=TRUE))

#Defin relevant working directories
data_dir<-"//storage.research.sesync.org/palmer-group-data/choptank/Nate/Storage_Capacity/"
scratch_dir<-"C:\\ScratchWorkspace\\"
wbt_dir<-"C:/WBT/whitebox_tools"

#Download packages 
library(tidyverse)
library(raster)
library(sf)
library(parallel)
library(devtools)

#Download package from GIT
install_github("FloodHydrology/InundationHydRology")
library(InundationHydRology)

#Download relevant data 
dem_jl<-raster(paste0(data_dir,"II_work/jl_dem"))
dem_jr<-raster(paste0(data_dir,"II_work/jr_dem"))
pnts<-readOGR(paste0(data_dir,"II_Work/."),"Wetland_Locations", verbose = F)
burn<-readOGR(paste0(data_dir,"II_Work/."),"wetland_burn")

#Add Unique ID to points
pnts$WetID<-seq(1,length(pnts))
burn$WetID<-sp::over(burn, pnts)$WetID
