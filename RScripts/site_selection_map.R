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
#2.0 Prep data==================================================================
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#crop areas to properties~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
mask<-st_buffer(site_bndry, 250)
subsheds<-subsheds[mask,]
wetlands<-giws[mask,]

#Create flowpath raster~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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
seeds<-wetlands %>% 
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

#Convert to bianary raster
flowpaths<-flowpaths %>% st_as_stars(.) %>% st_as_sf(., merge = TRUE)

#Checkout the results
flowpaths %>% st_transform(., crs=4326) %>% leaflet(.) %>% addPolygons(.)
  
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#3.0 Plot data==================================================================
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Fix output items for display
wetlands<- wetlands %>% 
  mutate(wetland_hsc_cm = if_else(is.na(wetland_hsc_cm), 
                                  merged_hsc_cm, 
                                  wetland_hsc_cm))


#Convert to coorindate systems [from planar]
site_bndry<-site_bndry %>% st_transform(.,crs=4326)
wetlands <- wetlands %>% st_transform(., crs = 4326)
connections <- flowpaths %>% st_transform(.,crs=4326)
wet_basins<-subsheds %>% st_transform(., crs=4326)

#color for leaflet map
pal <- colorBin("Blues", 
                domain = wetlands$watershed_hsc_cm, 
                bins = seq(0, base::max(wetlands$watershed_hsc_cm, na.rm=T), length.out = 9))

#Create lables
labels <- sprintf(
  "<strong>WetID = %s</strong>
   <br/>Watershed Area = %g ha
   <br/>Network Order = %g
   <br/>P/A Ratio = %g
   <br/>HSC [Wetland] = %g cm
   <br/>HSC [Watershed] = %g cm
   <br/>HAND = %g m
   <br/>HANS = %g m",
   wetlands$WetID, wetlands$watershed_area_m2/10000, wetlands$wet_order, wetlands$p_a_ratio, 
   wetlands$wetland_hsc_cm, wetlands$watershed_hsc_cm, wetlands$hand_m, wetlands$hans_m
) %>% lapply(htmltools::HTML)

#leaflet
m<-leaflet(wetlands) %>% 
  #Add Basemaps
  addProviderTiles("Esri.WorldImagery", group = "ESRI") %>% 
  addTiles(group = "OSM") %>%
  #Add properties
  addPolygons(data=site_bndry, weight =2, col="red", fill=F, group = "TNC Property") %>% 
  #Add subsheds
  addPolygons(data=wet_basins, weight =2, col="grey", group = "Subsheds") %>% 
  #Add flowpath data
  addPolygons(data=connections, weight=1, col="blue", group = "SW Connections") %>% 
  #Add wetlands
  addPolygons(
    fillColor = ~pal(watershed_hsc_cm), 
    weight = 2,
    opacity = 1,
    color = "grey90",
    dashArray = "3",
    fillOpacity = 0.7, 
    highlight = highlightOptions(
      weight = 5,
      color = "#666",
      dashArray = "",
      fillOpacity = 0.7,
      bringToFront = TRUE), 
    label = labels,
    labelOptions = labelOptions(
      style = list("font-weight" = "normal", padding = "3px 8px"),
      textsize = "15px",
      direction = "auto"), 
    group = "Wetlands") %>% 
  #Add Legend
  addLegend(pal = pal, values = ~watershed_hsc_cm, opacity = 0.7, title = NULL, position = "bottomright") %>% 
  #Add Layer Control Options
  addLayersControl(baseGroups = c("Esri", "OSM"), 
                   overlayGroups = c("Wetlands",
                                     "Subsheds",
                                     "SW Connections", 
                                     "TNC Property"
                                     ))

#Save map as html document
saveWidget(m, file="site_selection.html")
  
