#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Title: AIK Analysis
#Coder: C. Nathan Jones (cnjones7@ua.edu)
#Date: 8/30/2019
#Purpose: Estimate hydrogeomorphic features relevant to AIK's analysis
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#1.0 Setup workspace============================================================
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Clear memory
rm(list=ls(all=TRUE))

#Defin relevant working directories
data_dir<-"//nfs/palmer-group-data/choptank/Nate/Spatial_Analysis/"

#Download packages 
library(leaflet)
library(sf)
library(raster)
library(tidyverse)

#read data from spatial analysis
giws<-st_read(paste0(data_dir, "III_Products/giws.shp"))

#read anna's emperical data and join to spatial analysis
df<-read_csv(paste0(data_dir, "I_Data/q2_locations.csv")) %>% 
  st_as_sf(.,  coords = c("X", "Y"), crs = 4326) %>% 
  st_transform(., p) %>% 
  st_join(., giws)

#write output
#write_csv(df, paste0(data_dir, "III_Products/AIK_Metrics.csv"))