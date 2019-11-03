#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Title: Site Selection Map
#Coder: C. Nathan Jones (cnjones7@ua.edu)
#Date: 8/30/2019
#Purpose: Create Leaflet Map to assist with site selection
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

#For now, just do this
load('backup.RData')

#Define master projection [note, we have to reproject to coordinates later]
p<-"+proj=utm +zone=17 +ellps=GRS80 +datum=NAD83 +units=m +no_defs"

#Downlaod releavnt relavant layers
jr_shp<-st_read(paste0(data_dir,"/I_Data/JR_boundary.shp")) %>% st_transform(., p) %>% st_zm(.)
jl_shp<-st_read(paste0(data_dir,"/I_Data/JL_boundary.shp")) %>% st_transform(., p) %>% st_zm(.)
site_bndry<-st_union(jr_shp, jl_shp)

# #downlaod dem
dem<-raster(paste0(data_dir, "III_Products/dem_filter.tif"), overwrite=T)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#2.0 Filter sites based on group conversations==================================
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#crop areas to properties
mask<-st_buffer(site_bndry, 250)
subsheds<-subsheds[mask,]
giws<-giws[mask,]

#Select catchments with wetland orders between 3 and 8
sites<-giws %>% st_drop_geometry() %>% filter(wet_order<=6) %>% select(WetID) %>% as_vector(.)
watersheds<-watersheds %>% filter(WetID %in% sites)

#Limit giws to only "root" wetlands
giws<-giws %>% filter(WetID == root_giw)

#Estiamte catchment characteristics~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Create function to estimate catchment charactersitics
fun<-function(n){
  #Isolate catchment and giw
  watershed<-watersheds[n,]
  giw<-giws %>% filter(WetID == watershed$WetID)
  
  #estimate watershed storage capacity
  hsc_watershed_cm <- giw$watershed_hsc_cm
  
  #Estiamte wetland order
  wet_order<- giw$wet_order
  
  #esitmate number of wetlands
  wetlands<-giws[watershed,]
  n_wetlands<-  nrow(wetlands)
  
  #estimate gradient in HAND
  hand_range_m<-base::max(wetlands$hand_m)-base::min(wetlands$hand_m)
  
  #Create output
  output<-tibble(
    WetID = watershed$WetID,
    hsc_watershed_cm, 
    wet_order, 
    n_wetlands, 
    hand_range_m
  )
  
  #Export output
  output
}

#Apply function to existing catchments
output<-lapply(seq(1, nrow(watersheds)), fun) %>% bind_rows(.)

#Combine output with catchment shapes
watersheds<-left_join(watersheds, output)

#More filtering~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
watersheds <- watersheds %>% 
  filter(n_wetlands<6) %>% 
  filter(n_wetlands>2) %>% 
  filter(wet_order<=n_wetlands)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#3.0 Create flowpaths===================================================================
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#breach
breach_depressions(dem = paste0(data_dir, "III_Products/dem_filter.tif"), 
                   output = paste0(data_dir, "II_Work/breach.tif"), 
                   flat_increment = 0.01, 
                   fill_pits = T)

#FDR/FAC
d8_pointer(dem = paste0(data_dir, "II_Work/breach.tif"), 
           output = paste0(data_dir, "II_Work/fdr.tif"))
d8_flow_accumulation(dem = paste0(data_dir, "II_Work/breach.tif"), 
                     output = paste0(data_dir, "II_Work/fac.tif"))

#Define fac in R enviornment
fac<-raster(paste0(data_dir, "II_Work/fac.tif"))

#Define seed points
seeds<-giws %>% 
  #Filter to root giws only
  filter(is.na(merge_to)) %>%
  #Select unique ID for each polygon
  select(WetID) %>% 
  #Convert to points
  st_cast(., "MULTIPOLYGON") %>% 
  st_cast(., "MULTILINESTRING") %>% 
  st_cast(., "MULTIPOINT") %>% 
  st_cast(., "POINT") %>% 
  #Extract FAC raster value at each point
  mutate(fac = raster::extract(fac, .)) %>%
  #Select max fac point for each wetland
  group_by(WetID) %>%
  filter(fac == base::max(fac, na.rm=T))

#write seeds to workspace
st_write(seeds, paste0(data_dir, "II_Work/seeds.shp"), delete_layer = TRUE)

#Execute WBT flowpaths function
trace_downslope_flowpaths(seed_pts = paste0(data_dir,"II_Work/seeds.shp"), 
                          d8_pntr  = paste0(data_dir, "II_Work/fdr.tif"), 
                          output   = paste0(data_dir, "II_Work/flowpaths.tif"))

#Read flowpath raster into R environment
flowpaths<-raster(paste0(data_dir, "II_Work/flowpaths.tif"))

#Convert to bianary raster
flowpaths<-flowpaths*0+1

#Substract wetland areas
giws_grd<-(fasterize(giws, dem)*0)
giws_grd[is.na(giws_grd)]<-1
giws_grd[giws_grd==0]<-NA
flowpaths<-flowpaths*giws_grd

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#3.0 Create Map=================================================================
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Wetlands to trim:
bad<-c(903, 788, 493, 465, 921, 975, 919, 1124, 1152, 1029, 923, 718, 963, 965, 973, 1038, 979, 876, 800, 737, 689, 523, 509, 408, 1176, 978, 1001, 935, 1184, 1174, 1056, 727,  695, 504, 539)
watersheds<-watersheds %>% filter(!(WetID %in% bad))


#Project Layers
watersheds<-watersheds %>% st_transform(., crs=4326) 
giws<-giws %>% st_transform(., crs=4326) 
site_bndry<-site_bndry %>% st_transform(., crs=4326) 
flowpaths<-projectRaster(flowpaths, crs='+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs')

#Create lables
labels <- sprintf(
  "<strong>CatchmentID = %s</strong>
   <br/>Number of Wetlands = %g 
   <br/>Network Order = %g
   <br/>HSC [Watershed] = %g cm
   <br/>HAND Range = %g m",
  watersheds$WetID, watersheds$n_wetlands, watersheds$wet_order, 
  watersheds$hsc_watershed_cm, watersheds$hand_range_m
) %>% lapply(htmltools::HTML)

#leaflet
m<-leaflet() %>% 
  #Add Basemaps
  addProviderTiles("Esri.WorldImagery", group = "ESRI") %>% 
  addTiles(group = "OSM") %>%
  #Add properties
  addPolygons(data=site_bndry, 
              weight =2, 
              col="red", 
              fill=F, 
              group = "TNC Property") %>% 
  #Add flowpath data
  addRasterImage(x=flowpaths,
                 col="blue",
                 group = "SW Connections") %>%
  #Add wetlands
  addPolygons(
    data = giws, 
    weight = 1,
    opacity = 1,
    color = "blue") %>% 
  #Add watersheds
  addPolygons(
    data = watersheds, 
    weight = 2,
    color = "orange",
    fill=T,
    dashArray = "3",
    highlight = highlightOptions(
      weight = 5,
      color = "#666",
      dashArray = "",
      fillOpacity = 0,
      bringToFront = TRUE), 
    label = labels,
    labelOptions = labelOptions(
      style = list("font-weight" = "normal", padding = "3px 8px"),
      textsize = "15px",direction = "auto"), 
    group = "Watersheds") %>% 
  #Add Layer Control Options
  addLayersControl(baseGroups = c("Esri", "OSM"), 
                   overlayGroups = c("Watersheds",
                                     "SW Connections", 
                                     "TNC Property"
                   ))

#Export widget
saveWidget(m, file="site_selection_update.html")




