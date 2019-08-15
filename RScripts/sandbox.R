#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Title: Sandbox
#Coder: C. Nathan Jones (cnjones7@ua.edu)
#Date: 8/14/2019
#Purpose: Examine 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#1.0 Setup workspace------------------------------------------------------------
#Clear memory
rm(list=ls(all=TRUE))

#Defin relevant working directories
data_dir<-"//nfs/palmer-group-data/choptank/Nate/Spatial_Analysis/"

#Download packages 
library(gstat)
library(stars)
library(fasterize)
library(whitebox)
library(sf)
library(raster)
library(tidyverse)
library(parallel)

#Load custom R scripts
funs<-list.files("functions/")
for(i in 1:length(funs)){source(paste0('functions/',funs[i]))}
remove(funs)

#Define projection
p<-"+proj=utm +zone=17 +ellps=GRS80 +datum=NAD83 +units=m +no_defs"

#Download relevant data 
dem<-raster(paste0(data_dir,"I_Data/dem_jr"))
dem[dem<0]<-NA
dem<-projectRaster(dem, crs=p)  

#2.0 Filter the DEM-------------------------------------------------------------
#2.1 Fill Single Cell Pitts
#Export DEM to workspace
writeRaster(dem, paste0(data_dir,"II_Work/dem_test.tif"), overwrite=T)

#Fill single cell pits
fill_single_cell_pits(dem=paste0(data_dir,"II_Work/dem_test.tif"), 
                      output = paste0(data_dir,"II_Work/dem_fill.tif"))

#Read raster back into workspace
dem_fill<-raster(paste0(data_dir,"II_Work/dem_fill.tif"))

#2.2 Filter DEM
#Apply simple gausian filter to smooth random errors from DEM
dem_filter<- focal(dem_fill, w=focalWeight(dem, 3, "Gauss"))
crs(dem_filter)<-p

#3.0 Define "Root" Depressions--------------------------------------------------
#3.1 Use the Stochastic Deprresion Tool to identify deprresions
#Export fitlered DEM to workspace
writeRaster(dem_filter, paste0(data_dir,"II_Work/dem_filter.tif"), overwrite=T)

#Apply stochastic depressin analysis tool
stochastic_depression_analysis(dem = paste0(data_dir,"II_Work/dem_filter.tif"), 
                               output = paste0(data_dir,"II_Work/giws.tif"), 
                               rmse = 0.18, 
                               range = 10, 
                               iterations = 1000)

#3.2 Define depression based on threshold of occurence in stochastic procedure
#Reclass raster (any depression that was delineated less than 80% of the time is out!)
reclass(input = paste0(data_dir,"II_Work/giws.tif"), 
        output = paste0(data_dir,"II_Work/reclass.tif"), 
        reclass_vals = "'0;0;0.8'")
reclass(input = paste0(data_dir,"II_Work/reclass.tif"), 
        output = paste0(data_dir,"II_Work/reclass.tif"), 
        reclass_vals = "'1;0.8;1'")

#Convert 0 to NA
r<-raster(paste0(data_dir,"II_Work/reclass.tif")) 
r[r==0]<-NA
crs(r)<-p
writeRaster(r, paste0(data_dir, "II_Work/reclass_na.tif"), overwrite=T)
remove(r)

#Identify "clumps" of depressions
clump(input = paste0(data_dir, "II_Work/reclass_na.tif"), 
      output = paste0(data_dir, "II_Work/group.tif"), 
      zero_back = T)

#Export as polygon
giws<-raster(paste0(data_dir, "II_Work/group.tif")) 
giws[giws==1]<-NA
giws<- giws %>% st_as_stars(.) %>% st_as_sf(., merge = TRUE)

#Write polygon shapes to workspace
st_write(giws, paste0(data_dir, "II_Work/giws.shp"), delete_layer=TRUE)

#3.3 Filter depressions
#Filter by area and P:A Ratio
polygon_perimeter(input=paste0(data_dir, "II_Work/giws.shp"))
polygon_area(input=paste0(data_dir, "II_Work/giws.shp"))
giws<-st_read(paste0(data_dir, "II_Work/giws.shp"))
giws<-giws %>%
  #Remove small depressions
  filter(AREA>250) %>%
  #Remove oddly large depressions
  filter(AREA<1e5) %>%
  #Remove ditched depressions
  mutate(p_a_ratio = AREA/PERIMETER) %>%
  filter(p_a_ratio>2)

#3.4 Plot for funzies
plot(dem_filter)
plot(st_geometry(giws), add=T, col="blue")

#4.0 Define "leaf" or "child" depressions---------------------------------------
#For testing identify dragonfly bay complex
giw<-giws %>% filter(AREA==base::max(AREA))
n_slices<-10

#Create function to return conditional raster
con<-function(condition, trueValue, falseValue){
  return(condition * trueValue + (!condition)*falseValue)
}

#Crop DEM
temp<-crop(dem_filter, giw)
temp<-mask(temp, giw)

#Create Minimum Raster
temp_min<-temp*0+minValue(temp)
temp_min@crs<-temp@crs

#Estimate maximum depth
z_max<-cellStats(temp-temp_min, base::max)

#Create shape of maximum extent
max_inun<-con(temp>(temp_min+z_max), 0, 1) %>% st_as_stars() %>% st_as_sf(., merge = TRUE)

#Create levels tibble
nodes<-max_inun %>%
  mutate(node=1, 
         z = z_max,
         merge_to = NA, 
         spawned = 0) %>%
  select(node, z, merge_to, spawned)

#Drain the complex and identify spawns
for(i in 2:n_slices){
  #Define inundation depth
  z<-z_max-(z_max/n_slices)*i
  
  #define inundated areas
  inun<-con(temp>(temp_min+z),0,1)
  inun<-raster::clump(inun) %>% st_as_stars() %>% st_as_sf(., merge = TRUE)
  
  #Determine if any new nodes spawned 
  for(j in 1:nrow(nodes)){
    #Identify temp node
    node_temp<-nodes[j,]
    
    #Identify inundation within node
    node_inun<-inun[node_temp,]
    
    #If >2 shpae and the shape hasn't already spawned, then add spawn to node list
    if(nrow(node_inun)>1 & node_temp$spawned==0){
      #Add info to inundation shape
      node_inun<- node_inun %>%
        mutate(node = node_temp$node, 
               z = z, 
               merge_to = node_temp$node, 
               spawned = 0) %>%
        select(node, z, merge_to, spawned)
      
      #Add node id
      for(i in 1:nrow(node_inun)){
        node_inun$node[i]<-base::max(nodes$node)+i
      }
      
      #Add spawn indicator
      nodes$spawned[j]<-1
      
      #Merge with nodes
      nodes<-rbind(nodes, node_inun)
    }
  }
}

#Plot 
nodes %>% filter(!node %in% merge_to) %>% select(node) %>% plot()

#Next steps
# (1) Add filter for size
# (2) Create functin for individual depression



