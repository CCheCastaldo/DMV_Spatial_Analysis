#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Title: Site Selection Map
#Coder: C. Nathan Jones (cnjones7@ua.edu)
#Date: 3/5/2020
#Purpose: Create Leaflet Map to Display Well Locations
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#1.0 Setup workspace============================================================
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Clear memory
rm(list=ls(all=TRUE))

#Defin relevant working directories
data_dir<-"//nfs/palmer-group-data/choptank/Nate/Spatial_Analysis/"

#Download packages 
library(htmlwidgets)
library(leaflet)
library(whitebox)
library(stars)
library(fasterize)
library(sf)
library(raster)
library(tidyverse)

#Define master projection [note, we have to reproject to coordinates later]
p<-"+proj=utm +zone=17 +ellps=GRS80 +datum=NAD83 +units=m +no_defs"

#Downlaod boundary layers
jr_shp<-st_read(paste0(data_dir,"/I_Data/JR_boundary.shp")) %>% st_transform(., p) %>% st_zm(.)
jl_shp<-st_read(paste0(data_dir,"/I_Data/JL_boundary.shp")) %>% st_transform(., p) %>% st_zm(.)
site_bndry<-st_union(jr_shp, jl_shp)

#downlaod dem
dem<-raster(paste0(data_dir, "III_Products/dem_burn.tif"), overwrite=T)

#downlaod well lections
wells<-read_csv(paste0(data_dir,"I_Data/well_locations.csv")) %>% 
  st_as_sf(coords = c("POINT_X", "POINT_Y"), crs = 2248) %>% 
  st_transform(crs=p)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#2.0 Plot data==================================================================
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Fix output items for display
wetlands<- wetlands %>% 
  mutate(wetland_hsc_cm = if_else(is.na(wetland_hsc_cm), 
                                  merged_hsc_cm, 
                                  wetland_hsc_cm))


#Convert to coorindate systems [from planar]
site_bndry<-site_bndry %>% st_transform(.,crs=4326)
wells <- wells %>% st_transform(., crs = 4326)

#Create lables
labels <- sprintf("<strong>SiteID = %s</strong>",wells$Site_ID) %>% lapply(htmltools::HTML)

#leaflet
m<-leaflet(wells) %>% 
  #Add Basemaps
  addProviderTiles("Esri.WorldImagery", group = "ESRI") %>% 
  addTiles(group = "OSM") %>%
  #Add properties
  addPolygons(data=site_bndry, weight =2, col="red", fill=F, group = "TNC Property") %>% 
  #Add subsheds
  addCircleMarkers(data=wells, 
    weight =2, 
    col="orange",  
    label = labels,
    labelOptions = labelOptions(
      style = list("font-weight" = "normal", padding = "3px 8px"),
      textsize = "15px",
      direction = "auto"), 
    group = "Palmer Lab Wells") %>% 
  #Add Layer Control Options
  addLayersControl(baseGroups = c("Esri", "OSM"), 
                   overlayGroups = c("TNC Property",
                                     "Palmer Lab Wells"))

#Save map as html document
saveWidget(m, file="palmer_lab_well_locations.html")
  
m
