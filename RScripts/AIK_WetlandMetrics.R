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
#3) This anlaysis also relies on the InundationHydrology package, which has to be 
#     from directly from github. 

#To-do list:
#0) Update to latest WBT exe
#1) Use 2007 DEM to be consistent with KH
#2) Check with Anna to make sure Burn areas make sense
#3) Convert from sp hell to sf

#1.0 Setup workspace------------------------------------------------------------
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
pnts<-st_read(paste0(data_dir,"II_Work/AIK_Wetland_Locations.shp"))
burn<-st_read(paste0(data_dir,"II_Work/wetland_burn.shp"))

#Add Unique ID to points
pnts$WetID<-seq(1,nrow(pnts)) 
burn<-st_join(burn, pnts) %>% dplyr::select(WetID) %>% filter(!is.na(WetID))

#for now convert to sp (I'll need to update code to use sf at a later date)
pnts<-as_Spatial(pnts)
burn<-as_Spatial(burn)
library(rgeos)
library(rgdal)

#2.0 Burn Wetlands into DEM-----------------------------------------------------
#Process Jackson Ln property
mask_jl<-extent(dem_jl)
burn_jl<-crop(burn, mask_jl)
dem_burn_jl<-GIW_burn(dem_jl, burn_jl[1,])
for(i in 2:length(burn_jl)){
  dem_burn_jl<-GIW_burn(dem_burn_jl, burn_jl[i,])
}

#Process Jones Rd property
mask_jr<-extent(dem_jr)
burn_jr<-crop(burn, mask_jr)
dem_burn_jr<-GIW_burn(dem_jr, burn_jr[1,])
for(i in 2:length(burn_jr)){
  dem_burn_jr<-GIW_burn(dem_burn_jr, burn_jr[i,])
}

#3.0 Delineate Wetlands---------------------------------------------------------
#3.1 Identify deperessions------------------------------------------------------
#Identify depressions in DEM
giws_jl<-GIW_identification(
  dem=dem_burn_jl,                    #DEM in question
  min_size=100,                       #Minimum depression size (in map units)
  workspace="C:\\ScratchWorkspace\\", #Scratch Workspace
  wbt_path="C:/WBT/whitebox_tools")   #WBG toolbox location

giws_jr<-GIW_identification(
  dem=dem_burn_jr,                    #DEM in question
  min_size=100,                       #Minimum depression size (in map units)
  workspace="C:\\ScratchWorkspace\\", #Scratch Workspace
  wbt_path="C:/WBT/whitebox_tools")   #WBG toolbox location

#3.2 Wetland Subshed Delineation------------------------------------------------
#Create rapper function for subshed delineation
fun<-function(WetID){
  #Define point of interest
  pnt<-pnts[pnts$WetID==WetID,]
  
  #Define dem and depressions of interest
  if(gIntersects(pnt, burn_jr)==T){
    dem=dem_burn_jr
    depressions=giws_jr
  }else{
    dem=dem_burn_jl
    depressions=giws_jl
  }
  
  #Execute subshed delineation function
  subshed<-GIW_subshed_delineation(
    dem = dem,
    depressions=depressions,
    wetland=pnt)
  
  #Write raster to folder
  writeRaster(subshed, 
              paste0(data_dir, "II_work/subshed_",pnt$WetID,".asc"),
              overwrite=T)
}

#Run function
lapply(seq(1,nrow(burn)), fun)
