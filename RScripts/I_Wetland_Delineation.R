#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Title: Spatial Analysis
#Coder: C. Nathan Jones (cnjones7@ua.edu)
#Date: 11/3/2019
#Purpose: Examine hydrogeomorphic features accross Delmarva Bay wetland sites
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#1.0 Setup workspace============================================================
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Clear memory
rm(list=ls(all=TRUE))

#Defin relevant working directories
data_dir<-"//nfs/palmer-group-data/choptank/Nate/Spatial_Analysis/"

#Download packages 
library(gstat)
library(stars)
library(fasterize)
library(whitebox)
library(igraph)
library(sf)
library(raster)
library(tidyverse)
library(parallel)

#Define projection
p<-"+proj=utm +zone=17 +ellps=GRS80 +datum=NAD83 +units=m +no_defs"

#Download relevant data 
dem<-raster(paste0(data_dir,"I_Data/2007_1m_DEM.tif"))
  dem[dem<0]<-NA
  dem<-projectRaster(dem, crs=p)  
streams<-st_read(paste0(data_dir, "I_Data/TotalStreamNetwork_UCRW.shp"))
  streams<-st_transform(streams, crs=p)
  streams<-st_crop(streams, dem)
burn<-st_read(paste0(data_dir, "I_Data/FN_GB_Wetland_Burn.shp"))
  burn<-st_transform(burn, crs=p)
  
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#2.0 Prep DEM for Analysis======================================================
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#2.1 Filter the DEM-------------------------------------------------------------
#Export DEM to workspace
writeRaster(dem, paste0(data_dir,"II_Work/dem.tif"), overwrite=T)

#Fill single cell pits
fill_single_cell_pits(dem=paste0(data_dir,"II_Work/dem.tif"), 
                      output = paste0(data_dir,"II_Work/dem_fill.tif"))

#Read raster back into workspace
dem_fill<-raster(paste0(data_dir,"II_Work/dem_fill.tif"))

#Apply simple gausian filter to smooth random errors from DEM
dem_filter<- focal(dem_fill, w=focalWeight(dem, 3, "Gauss"))
crs(dem_filter)<-p

#2.2 Burn Wetlands of interest into DEM-----------------------------------------
#Create DEM urn Function~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
GIW_burn<-function(
  dem,   #DEM in question 
  burn){ #Wetland Polygon
  
  #Create DEM mask
  mask<-rasterize(burn, dem, 1)
  dem_mask<-mask*dem
  
  #Create minimum raster
  dem_min<-cellStats(dem_mask, base::min, na.rm=T)
  dem_min<-dem_mask*0+dem_min
  dem_min[is.na(dem_min)]<-0
  
  #Replace masked location with min raster value
  dem_mask<-dem_mask*0
  dem_mask[is.na(dem_mask)]<-1
  dem<-dem*dem_mask+dem_min
  
  #Export DEM
  dem
}

#Apply function to wetlands of interest
dem_burn<-GIW_burn(dem=dem_filter, burn=burn[1,])
dem_burn<-GIW_burn(dem=dem_filter, burn=burn[2,])

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#3.0 Define "Root" Depressions==================================================
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#3.1 Use the Stochastic Deprresion Tool to identify deprresions-----------------
#Export fitlered DEM to workspace
writeRaster(dem_burn, paste0(data_dir,"II_Work/dem_burn.tif"), overwrite=T)

#Apply stochastic depressin analysis tool
set.seed(100)
stochastic_depression_analysis(dem = paste0(data_dir,"II_Work/dem_burn.tif"), 
                               output = paste0(data_dir,"II_Work/giws.tif"), 
                               rmse = 0.18, 
                               range = 10, 
                               iterations = 100)

#3.2 Define depression based on threshold of occurence in stochastic procedure----
#Reclass raster (any depression that was delineated less than 80% of the time is out!)
reclass(input = paste0(data_dir,"II_Work/giws.tif"), 
        output = paste0(data_dir,"II_Work/reclass.tif"), 
        reclass_vals = "'0;0;0.8'")
reclass(input = paste0(data_dir,"II_Work/reclass.tif"), 
        output = paste0(data_dir,"II_Work/reclass.tif"), 
        reclass_vals = "'1;0.8;1'")

#Convert 0 to NA
giws<-raster(paste0(data_dir,"II_Work/reclass.tif")) 
giws<-raster::clump(giws)

#Export as polygon
giws[giws==1]<-NA
giws<- giws %>% st_as_stars(.) %>% st_as_sf(., merge = TRUE)

#Write polygon shapes to workspace
st_write(giws, paste0(data_dir, "II_Work/giws.shp"), delete_layer=TRUE)

#3.3 Filter depressions---------------------------------------------------------
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

#3.4 Add Unique Identifier to each GIW------------------------------------------ 
giws$WetID<-seq(1, nrow(giws))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#4.0 Define merging patterns within depressions=================================
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Reference to Wu et al [2019] code: 
#https://github.com/giswqs/lidar/blob/master/lidar/slicing.py

#4.1 Create function to execute level set method individual basins--------------
fun<-function(n, giws, dem_burn, max_size = 250, n_slices = 10){
  
  #Identify GIW
  giw<-giws[n,]
  
  #Create function to return conditional raster
  con<-function(condition, trueValue, falseValue){
    return(condition * trueValue + (!condition)*falseValue)
  }
  
  #Crop DEM
  temp<-crop(dem_burn, giw)
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
    mutate(WetID=giw$WetID, 
           z = z_max,
           merge_to = NA, 
           spawned = 0) %>%
    select(WetID, z, merge_to, spawned)
  
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
      
      #Filter out small wetlands
      node_inun<-node_inun %>%
        mutate(size=as.numeric(st_area(node_inun))) %>%
        filter(size>max_size)
      
      #If >2 shpae and the shape hasn't already spawned, then add spawn to node list
      if(nrow(node_inun)>1 & node_temp$spawned==0){
        #Add info to inundation shape
        node_inun<- node_inun %>%
          mutate(WetID = node_temp$WetID, 
                 z = z, 
                 merge_to = node_temp$WetID, 
                 spawned = 0) %>%
          select(WetID, z, merge_to, spawned)
        
        #Add node id
        for(k in 1:nrow(node_inun)){
          node_inun$WetID[k]<-paste0(giw$WetID,"-",i,"-",j,"-",k)
        }
        
        #Add spawn indicator
        nodes$spawned[j]<-1
        
        #Merge with nodes
        nodes<-rbind(nodes, node_inun)
      }
    }
  }
  
  #Add infomration about root wetland
  nodes$root_giw<-giw$WetID
  
  #Export nodes
  nodes
}

#4.2 Apply function-------------------------------------------------------------
outside_fun<-function(x){
  tryCatch(fun(x, giws, dem_burn, max_size = 250, n_slices = 10), 
           error = function(e) NA)
}

#apply function (~ minutes on SESYNC server)
t0<-Sys.time()
giws<-mclapply(X=seq(1,nrow(giws)), FUN=outside_fun, mc.cores=detectCores())
giws<-do.call(rbind, giws)
tf<-Sys.time()
tf-t0

#4.3 Convert Unique GIW IDs to numeric format-----------------------------------
#Create "change" tibble
index<-giws %>% 
  #Highlight ID of branch and leaf wetlands (i.e., ones with "-")
  filter(str_detect(WetID, "-")) %>%
  #create new list of new IDs
  mutate(NewID = seq(base::max(as.numeric(paste(giws$WetID)), na.rm=T)+1,  #Max WetID + 1 
                     base::max(as.numeric(paste(giws$WetID)), na.rm=T)+nrow(.))) %>%  #Max WetID + number of spawned wetlands
  #Make NewID a character to mach WetID [for now!]
  mutate(NewID=paste(NewID)) %>%
  #Select old and new IDs
  as_tibble() %>%
  select(WetID, NewID)

#Update values in WetID collumn
giws<-giws %>% 
  #left join index
  left_join(., index) %>% 
  #conditionally chnage WetID value
  mutate(WetID = if_else(is.na(NewID), WetID, NewID)) %>%
  #Convert WetID to numeric
  mutate(WetID = as.numeric(paste(WetID))) %>%
  #Remove NewID Collumn
  select(-NewID)

#Update merge_to collumn
index<-index %>% rename(merge_to = WetID)
giws<-giws %>% 
  #left join index
  left_join(., index) %>% 
  #conditionally chnage WetID value
  mutate(merge_to = if_else(is.na(NewID), merge_to, NewID)) %>%
  #Convert WetID to numeric
  mutate(merge_to = as.numeric(paste(merge_to))) %>%
  #Remove NewID Collumn
  select(-NewID)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#5.0 Define connectivity between wetlands=======================================
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#5.1 Create flow accumulation and flow direction rasters------------------------
#Export raster to workspace
writeRaster(dem_burn, paste0(data_dir, "II_Work/dem.tif"), overwrite=T)

#Execute breach tool
breach_depressions(dem = paste0(data_dir,"II_Work/dem.tif"), 
                   output = paste0(data_dir,"II_Work/dem_breach.tif"), 
                   fill_pits = TRUE)

#Execute WBT flow accumulation raster
d8_flow_accumulation(dem = paste0(data_dir,"II_Work/dem_breach.tif"), 
                     output = paste0(data_dir,"II_Work/fac.tif"))

d8_pointer(dem = paste0(data_dir, "II_Work/dem_breach.tif"), 
           output = paste0(data_dir, "II_Work/fdr.tif"))

#Pull fac back into R environment
fac<-raster(paste0(data_dir, "II_Work/fac.tif"))
fdr<-raster(paste0(data_dir, "II_Work/fdr.tif"))

#5.2 Create functin to identify "down gradient" wetland-------------------------
fun<-function(n,
              giws,
              fac,
              fdr){
  
  #Identify giw of interest 
  giw<- giws[n,]
  
  #If a leaf or branch wetland, skip!
  if(is.na(giw$merge_to)){
    
    #Create temporary workspace
    scratch_dir<-paste0(tempfile(),"/")
    dir.create(scratch_dir)
    
    #Find max fac point within giw
    fac_giw<-crop(fac, as_Spatial(giw))
    fac_giw<-mask(fac_giw, as_Spatial(giw))
    
    #create pour point
    pp<-rasterToPoints(fac_giw) %>% 
      #convert to tibble
      as_tibble(.) %>%
      #filter to just max fac value
      filter(fac == base::max(fac, na.rm=T)) %>%
      #Select first row
      slice(1) %>%
      #Make pour point an sf shape
      st_as_sf(., coords = c("x", "y"), crs = st_crs(giw))
    
    #Export files to scratch dir
    writeRaster(fdr, paste0(scratch_dir, "fdr.tif"), overwrite=T)
    write_sf(pp, paste0(scratch_dir,"pp.shp"), delete_layer = T)
    
    #Execute WBT flowpaths functin
    trace_downslope_flowpaths(seed_pts = paste0(scratch_dir,"pp.shp"), 
                              d8_pntr  = paste0(scratch_dir, "fdr.tif"), 
                              output   = paste0(scratch_dir, "flowpath.tif"))
    
    #Read flowpath raster into R workspace
    flowpath_grd<-raster(paste0(scratch_dir, "flowpath.tif")) 
    
    #Convert flowpath raster to sf shape
    flowpath_shp <- flowpath_grd %>% st_as_stars() %>% st_as_sf(., merge = TRUE)
    
    #Identify downstream giws (root GIWS only)
    giws_ds<-giws[flowpath_shp,] %>% filter(is.na(merge_to)) %>% filter(WetID != giw$WetID)
    
    #Convert to points
    giw_bndry<- giws_ds %>% select(WetID) %>% st_segmentize(.,1) %>% st_cast(., "MULTILINESTRING") %>% st_cast(.,"MULTIPOINT") %>% st_cast(.,"POINT")
    
    #Estimate flow accumulation along flowpath
    fac_fp<-flowpath_grd*fac
    
    #Estimate downstream wetland [i.e., shape with minimum flowpath fac value]
    if(nrow(giws_ds)>0){
      output<-giw_bndry %>%
        #Extract flowpath fac data for downstream giw boundaries
        mutate(fac = raster::extract(fac_fp, .)) %>%
        #select minimum point 
        filter(fac == base::min(fac, na.rm=T)) %>% slice(1) %>%
        #Create output collumns 
        mutate(flow_to = WetID, 
               WetID = giw$WetID) %>% 
        #Convert to tibble
        as_tibble() %>%
        #Select collumns of interest
        select(WetID, flow_to) 
    }else{
      output<-tibble(WetID = giw$WetID, 
                     flow_to = NA)
    }
    
    #Delete temp file
    unlink(scratch_dir, recursive = T)
    
  }else{
    output<-tibble(WetID = giw$WetID, 
                   flow_to = NA)
  }
  
  #Export output
  output
}

#5.3 Apply function-------------------------------------------------------------
outside_fun<-function(x){
  tryCatch(fun(x, giws, fac, fdr), 
           error = function(e) NA)
}

#apply function (~ minutes on SESYNC server)
t0<-Sys.time()
output<-mclapply(X=seq(1,nrow(giws)), FUN=outside_fun, mc.cores=detectCores())
output<-do.call(rbind, output)
tf<-Sys.time()
tf-t0

#Merge output with giws tibble
giws<-left_join(giws, output)

#5.4 Create flowpath lines for plotting-----------------------------------------
#Define seed points
seeds<- giws %>% 
  #Filter to root giws only
  filter(is.na(merge_to)) %>%
  #Select unique ID for each polygon
  select(WetID) %>% 
  #Convert to points
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
                          output   = paste0(data_dir, "II_Work/flowpath.tif"))

#Read flowpath raster into R environment
flowpath<-raster(paste0(data_dir, "II_Work/flowpath.tif"))

#Convert to bianary raster
flowpath<-flowpath*0+1

#Substract wetland areas
giws_grd<-(fasterize(giws, dem)*0)
giws_grd[is.na(giws_grd)]<-1
giws_grd[giws_grd==0]<-NA
flowpath<-flowpath*giws_grd

#Convert to vector
flowpath<-flowpath %>% st_as_stars(.) %>% st_as_sf(., merge = TRUE)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#6.0 Watershed Delineation======================================================
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#6.1 Define pour points for all root and leaf wetlands--------------------------
#6.1.1 Create flowpath raster (yes, this is replicate from above)~~~~~~~~~~~~~~~
#Define seed points
seeds<- giws %>% st_cast(., "POINT") 

#write seeds to workspace
st_write(seeds, paste0(data_dir, "II_Work/seeds.shp"), delete_layer = TRUE)

#Execute WBT flowpaths function
trace_downslope_flowpaths(seed_pts = paste0(data_dir,"II_Work/seeds.shp"), 
                          d8_pntr  = paste0(data_dir, "II_Work/fdr.tif"), 
                          output   = paste0(data_dir, "II_Work/flowpath.tif"))

#Read flowpath raster into R environment
fp<-raster(paste0(data_dir, "II_Work/flowpath.tif"))*0+1

#6.1.2 Define pourpoint along flowpath for non-branch wetland~~~~~~~~~~~~~~~~~~~
#Add fac information to flowpath raster
fp<-fac*fp

#Create pp 
giws_grd<- giws %>%
  #Select non "branch" wetland shpaes
  filter(is.na(merge_to) | (merge_to>0 & spawned == 0)) %>%
  #Create raster where values correspond to WetID
  fasterize(., raster = dem, field = 'WetID')

#Estiamte pourpoint 
pp<-fp %>% 
  #Convert to sf points object
  rasterToPoints(.) %>%
  as_tibble(.) %>%
  st_as_sf(., coords = c('x','y'), crs=p) %>% 
  rename(fac=layer) %>%
  #Add WetID data
  mutate(WetID = raster::extract(giws_grd, .))%>%
  #Select the max fac per WetID
  na.omit() %>%
  group_by(WetID) %>%
  filter(fac == base::max(fac, na.rm=T)) %>%
  select(WetID)

#6.2 Delineate subsheds---------------------------------------------------------
#6.2.1 Subshed Delineation~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Write shpaes to workspace
st_write(pp, paste0(data_dir, "II_Work/pp.shp"), delete_layer=T)
writeRaster(dem_burn, paste0(data_dir, "II_Work/dem.tif"), overwrite=T)

#Conduct watershed analysis
watershed(d8_pntr = paste0(data_dir, "II_Work/fdr.tif"), 
          pour_pts = paste0(data_dir,"II_Work/pp.shp"), 
          output = paste0(data_dir,"II_Work/subshed.tif"))

#Read subsheds into R environment
subsheds<-raster(paste0(data_dir,"II_Work/subshed.tif"))

#Convert to polygon
subsheds<-subsheds %>% st_as_stars(.) %>% st_as_sf(., merge=T) 

#6.2.2 Add WetID to subsheds~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Create function to identify subhsed
fun<-function(n,
              subsheds,
              pp,
              fp){
  
  #Identify pour point of interest'
  p<-pp[n,]
  
  #Clip flowpath to pour point neighborhood
  flowpath_clip<-crop(fp, st_buffer(p,res(fp)[1]*5))
  flowpath_clip<-mask(flowpath_clip, st_buffer(p,res(fp)[1]*5))
  
  #Extract fac value at pp
  p_fac_value<-raster::extract(flowpath_clip, p)
  
  #Make raster cell at pp equal to NA
  flowpath_clip[flowpath_clip==p_fac_value]<-NA
  
  #Find "upstream" cell of pourpoint
  pp_new<-rasterToPoints(flowpath_clip) %>%
    #Convert to tibble for processing
    as_tibble(.) %>%
    #Remove "downstream points"
    filter(layer <= p_fac_value) %>%
    #Select max point
    filter(layer == base::max(layer)) %>%
    #Conver to sf point
    st_as_sf(.,coords = c("x", "y"), crs = st_crs(pp))
  
  #Identify overlapping subshed
  subshed<-subsheds[pp_new,]
  
  #Create output of subhsed id and WetID
  output<-tibble(
    WetID = p$WetID, 
    subshed = subshed$subshed
  )
  
  #Export output
  output
}

#Apply function
outside_fun<-function(x){
  tryCatch(fun(x, subsheds, pp, fp), 
           error = function(e) NA)
}

#apply function (~ minutes on SESYNC server)
t0<-Sys.time()
output<-mclapply(X=seq(1,nrow(pp)), FUN=outside_fun, mc.cores=detectCores())
output<-bind_rows(output)
tf<-Sys.time()
tf-t0

#Add WetID to subsheds
subsheds<-left_join(subsheds, output)

#6.2.3 Correct "merged" subsheds~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#We delineated leaf and root subsheds.  The leaf subsheds were deleted from the large 
#merged subshed, so we need to add them back in.

#Create list of merged subsheds
merge_id<- giws %>% as_tibble() %>% filter(root_giw != WetID) %>% select(root_giw) %>% unique() %>% as_vector() %>% as.numeric()

#create function to merge shapes (there's likely a better way to do this with purr)
fun<-function(n,
              merge_id,
              giws){
  #Identify subhseds to merge
  merged_subsheds<-
    #Convert giws to tibble
    giws %>% as_tibble() %>% 
    #Select wetlands in merged complex
    filter(root_giw == as.numeric(merge_id[n])) %>%
    #Reduce tibble to just the unique id
    select(WetID) %>%
    #left join with subshed
    left_join(.,subsheds) %>% 
    #Convert to a single polygon
    st_as_sf() %>% st_union() %>% st_sf() %>%
    #Add required collumns
    mutate(WetID = as.numeric(merge_id[n]), subshed = NA)
  
  #Export merged shape
  merged_subsheds
}

#Create update subshed shapes (i.e., apply the damn function)
outside_fun<-function(x){
  tryCatch(fun(x, merge_id, giws), 
           error = function(e) NA)
}
updated_subsheds<-lapply(seq(1,length(merge_id)), outside_fun)  %>% do.call(rbind, .)

#Replace "bad" subsheds with updated subsheds
subsheds<- subsheds %>% filter(!(WetID %in% merge_id)) %>% rbind(., updated_subsheds)

#6.3 Delineate watersheds-------------------------------------------------------
#Unnest subsheds using igraph approach. (Credit Goes to Kelly Hondula: https://github.com/ecodasXIII/coastalsheds/blob/master/02-get-huc12shed.Rmd#L105)
fun<-function(n, giws, subsheds){
  #Identify unique ID for wetland of interest (note, here we are choosing the root wetland's ID) 
  UID<-giws$root_giw[n] 
  
  #Create edgelist
  edgelist<-giws %>% 
    st_drop_geometry() %>% 
    select(WetID,flow_to) 
  
  #Create network object from edgelist
  network<-edgelist %>% graph_from_data_frame()
  
  #define upstream paths 
  paths_in<-all_simple_paths(network, from = paste(UID), mode = "in")
  upstream_subsheds<-sapply(paths_in, names) %>% unlist() %>% unique() %>% as.numeric()
  
  #Add focal wetlands subshed ID
  upstream_subsheds<-c(UID, upstream_subsheds)
  
  #Merge watershed shape
  watershed<-
    #Convert giws to tibble
    giws %>% as_tibble() %>% 
    #Select upstream wetlands
    filter(WetID %in% upstream_subsheds) %>%
    #Reduce tibble to just the unique id
    select(WetID) %>%
    #left join with subshed
    left_join(.,subsheds) %>% 
    #Convert to a single polygon
    st_as_sf() %>% st_union() %>% st_sf() %>%
    #Add unique identifier [note, use the WetID instead of root_ID here]
    mutate(WetID = giws$Wet[n])
  
  #Export watershed shape
  watershed
}

#Apply function
outside_fun<-function(x){
  tryCatch(fun(x, giws, subsheds), 
           error = function(e) NA)
}

#apply function (~ minutes on SESYNC server)
watersheds<-mclapply(X=seq(1,nrow(giws)), FUN=outside_fun, mc.cores=detectCores()) %>% do.call(rbind,.)

#6.4 Estimate area metrics------------------------------------------------------
#6.4.1 Estimate subshed area~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
subshed_area<-subsheds %>%
  #estimate area
  mutate(subshed_area_m2 = st_area(., by_element=T)) %>%
  #Get rid of units
  mutate(subshed_area_m2 = as.numeric(subshed_area_m2)) %>%
  #Convert to tibble
  as_tibble(.) %>%
  select(WetID, subshed_area_m2)

#Join to GIWS tibble
giws<-left_join(giws, subshed_area)

#6.4.2 Estimate watershed area area~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
watershed_area<-watersheds %>%
  #estimate area
  mutate(watershed_area_m2 = st_area(., by_element=T)) %>%
  #Get rid of units
  mutate(watershed_area_m2 = as.numeric(watershed_area_m2)) %>%
  #Convert to tibble
  as_tibble(.) %>%
  select(WetID, watershed_area_m2)

#Join to GIWS tibble
giws<-left_join(giws, watershed_area) 

#Filter out "slivers"
giws<-giws %>% 
  group_by(WetID) %>% 
  filter(subshed_area_m2 == base::max(subshed_area_m2)) %>% 
  ungroup(.)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#7.0 Estimtate Storage Capacity [Volume]========================================
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#7.1 Create function to estimate storage capacity-------------------------------
fun<-function(n, 
              giws, 
              dem){
  
  #Identify giw of interest
  giw<-giws[n,]
  
  #Convert to raster
  dem_giw<-crop(dem, giw)
  dem_giw<-mask(dem_giw, giw)
  
  #Create max raster
  max_giw<-dem_giw*0+cellStats(dem_giw, base::max)
  
  #Estimate volume between max_giw and dem_giw rasters
  depth_giw<-max_giw - dem_giw
  volume_m3<-cellStats(depth_giw, base::sum)*res(depth_giw)[1]*res(depth_giw)[2]
  
  #Export output
  tibble(WetID = giw$WetID, 
         volume_m3)
}


#7.2 Execute function-----------------------------------------------------------
outside_fun<-function(x){
  tryCatch(fun(x, giws, dem), 
           error = function(e) tibble(WetID=giws$WetID[n], 
                                      volume_m3=NA)
  )
}

#apply function (~ minutes on SESYNC server)
t0<-Sys.time()
output<-mclapply(X=seq(1,nrow(giws)), FUN=outside_fun, mc.cores=detectCores())
output<-do.call(rbind, output)
tf<-Sys.time()
tf-t0

#Merge output with giws tibble
giws<-left_join(giws, output)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#8.0 HSC Metrics================================================================
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#8.1 Individual Wetlands--------------------------------------------------------
#Estimate hsc for "leaf" wetlands
output<-giws %>% as_tibble() %>% 
  #Select lowest subunit of wetland (i.e.,a leaf)
  filter(spawned==0) %>%
  #Estimate Storage Capacity
  mutate(wetland_hsc_cm = volume_m3/subshed_area_m2*100) %>%
  #Select collumns for output
  select(WetID, wetland_hsc_cm)

#Join to GIWs tibble
giws<-left_join(giws, output)

#8.2 Root/Merged Wetlands-------------------------------------------------------
#Estimate hsc for merged wetlands [only max merge]
output<-giws %>% as_tibble() %>% 
  #Select lowest subunit of wetland (i.e.,a leaf)
  filter(spawned!=0 & WetID==root_giw) %>%
  #Estimate Storage Capacity
  mutate(merged_hsc_cm = volume_m3/subshed_area_m2*100) %>%
  #Select collumns for output
  select(WetID, merged_hsc_cm)

#Join to GIWs tibble
giws<-left_join(giws, output)

#8.3 Watersheds-----------------------------------------------------------------
#Creat function to find wetlands within a given wateshed
fun<-function(n, giws){
  
  #Select of interest 
  giw<-giws[n,]
  
  #Define its WetID [for later use]
  WetID<-giw$WetID
  
  #Identify root wetland
  #Is this a root wetlands?
  if(giw$WetID != giw$root_giw){
    #If not, redefine giw as its root wetland
    giw<-giws %>% filter(WetID == giw$root_giw)
  }
  
  #Identify upstream giws
  edgelist<-giws %>% 
    st_drop_geometry() %>% 
    select(WetID,flow_to) 
  network<-edgelist %>% graph_from_data_frame()
  paths_in<-all_simple_paths(network, from = paste(giw$WetID), mode = "in")
  upstream_giws<-sapply(paths_in, names) %>% unlist() %>% unique() %>% as.numeric()
  
  #add current giw to upstream giws tibble
  upstream_giws<-c(giw$WetID, upstream_giws)
  
  #Estiamte storage capacity in watershed
  watershed_hsc_cm<-giws %>%
    #Drop spatial attributes 
    st_drop_geometry() %>% 
    #filter to upstream giws
    filter(WetID %in% upstream_giws) %>% 
    #select volume
    select(volume_m3, watershed_area_m2) %>% 
    #sum
    summarise(watershed_hsc_cm = sum(volume_m3, na.rm=T)/giw$watershed_area_m2)*100 
  
  #Create Output
  output<-tibble(WetID, 
                 watershed_hsc_cm= watershed_hsc_cm$watershed_hsc_cm)
  output
}

#Create wrapper function w/ tryCatch for error handling
outside_fun<-function(x){
  tryCatch(fun(x, giws=giws), 
           error = function(e) tibble(WetID = giws$WetID[n],
                                      watershed_hsc_cm = NA))
}

#apply function 
output<-mclapply(X=seq(1,nrow(giws)), FUN=outside_fun, mc.cores=detectCores())
output<-bind_rows(output) %>% distinct(.)

#Join to GIWs tibble
giws<-left_join(giws, output)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#9.0 Estimate shape metrics=====================================================
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Write giw shape to working dir
giws %>% 
  select(WetID) %>% 
  st_write(., 
           paste0(data_dir, "II_Work/giws.shp"), 
           delete_layer = T)

#Estimate Long Axis wtih WBT
polygon_long_axis(input=paste0(data_dir, "II_Work/giws.shp"), 
                  output = paste0(data_dir, "II_Work/long.shp"))

#Read long axis shape into R env and
giws<-st_read(paste0(data_dir, "II_Work/long.shp")) %>% 
  #Estimate length
  mutate(a_axis_length_m = as.numeric(st_length(.))) %>% 
  #Drop geo and add to giws tibble
  st_drop_geometry() %>% left_join(giws, .)


#Estimate short Axis wtih WBT
polygon_short_axis(input=paste0(data_dir, "II_Work/giws.shp"), 
                   output = paste0(data_dir, "II_Work/short.shp"))

#Read short axis shape into R env and
giws<-st_read(paste0(data_dir, "II_Work/short.shp")) %>% 
  #Estimate length
  mutate(a_axis_length_m = as.numeric(st_length(.))) %>% 
  #Drop geo and add to giws tibble
  st_drop_geometry() %>% left_join(giws, .)


#Estimate shape parameters with WBT
polygon_perimeter(input=paste0(data_dir, "II_Work/giws.shp"))
polygon_area(input=paste0(data_dir, "II_Work/giws.shp"))
giws<-st_read(paste0(data_dir, "II_Work/giws.shp")) %>% 
  #Drop geo
  st_drop_geometry() %>% 
  #Rename WBT output
  rename(perimeter_m = PERIMETER, 
         area_m2 = AREA) %>% 
  #Estimate P:A Ratio
  mutate(p_a_ratio = perimeter_m/area_m2) %>% 
  #merge with GIW tibble
  left_join(giws, .)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#10.0 Estimate Height above nearest drainage (HAND) Metrics======================
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Convert stream network to raster
stream_grd<-streams %>% 
  #Create negative buffer
  st_buffer(., -1) %>% 
  #FASTERIZE!
  fasterize(., dem_burn)

#Write to workspace
writeRaster(stream_grd, paste0(data_dir, "II_Work/stream.tif"), overwrite = T)
writeRaster(dem_burn, paste0(data_dir,"II_Work/dem_burn.tif"), overwrite=T)

#Conduct breach analysis
breach_depressions(dem = paste0(data_dir,"II_Work/dem_burn.tif"), 
                   fill_pits = T, 
                   flat_increment = 0.001, 
                   output = paste0(data_dir,"II_Work/dem_breach.tif"))

#Execute elevation above stream tool
elevation_above_stream(
  dem= paste0(data_dir,"II_Work/dem_breach.tif"), 
  streams = paste0(data_dir, "II_Work/stream.tif"), 
  output = paste0(data_dir, "II_Work/hand.tif")
)

#Bring hand raster into R environment
hand_grd<-raster(paste0(data_dir, "II_Work/hand.tif"))

#Crop out unrealistic edge values
hand_grd[hand_grd>100]<-NA

#Crop out stream grid
stream_grd<-streams %>% fasterize(., dem_burn)
stream_grd<-stream_grd*0
stream_grd[is.na(stream_grd)]<-1
stream_grd[stream_grd==0]<-NA
hand_grd<-hand_grd*stream_grd

#Estimate hand for all giws
output<-mclapply(
  #input to function
  seq(1, nrow(giws)), 
  #function
  FUN = function(n){
    giw<-giws[n,]
    hand_temp<-crop(hand_grd, giw)
    hand_temp<-mask(hand_temp, giw)
    tibble(WetID = giw$WetID, 
           hand_m = cellStats(hand_temp, mean))
  }, 
  #cores
  mc.cores = detectCores()
) 

#Bind rows and merge with giws
giws<-output %>% bind_rows(.) %>% left_join(giws,.)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#11.0 Estimate Height Above Neighborhood Storage (HANS) metrics=================
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#11.1 Estiamte mean depth for each wetland--------------------------------------
fun<-function(n){
  #Isolate DEM
  giw<-giws[n,]
  
  #Crop DEM
  giw_dem<-crop(dem, giw)
  giw_dem<-mask(giw_dem, giw)
  
  #Estimate cell stats
  output = tibble( 
    WetID = giw$WetID,
    mean_elevation_m = cellStats(giw_dem, base::mean)
  )
}

#Apply function
output<-mclapply(X=seq(1, nrow(giws)), 
                 FUN = fun, 
                 mc.cores=detectCores()) %>% bind_rows()

#Join to giw tibble
giws<-left_join(giws, output)

#11.2 Create function to estimate HANS------------------------------------------
fun<-function(n, giws){
  
  #Isolate wetland of interest
  giw<-giws[n,] %>% st_centroid(.)
  
  #Creat mask with 1km^2 area 
  mask<-giw %>% st_buffer(., 564.1896) 
  
  #Select giws of interest
  giws_masked<-giws[mask,] %>% st_centroid(.) %>% filter(WetID!=giw$WetID)
  
  #Estimate distance between giws
  giws_masked$dist_m <- as.numeric(st_distance(giw, giws_masked, by_element = T))
  
  #Estimate weighted average wetland depth
  hans<-giws_masked %>% 
    #drop_geo
    st_drop_geometry(.) %>% 
    #Estimate inverse distance and elevation:distance ratio
    mutate(ele_dis = mean_elevation_m/dist_m, 
           inv_dis = 1/dist_m) %>% 
    #Sumarise IDW
    summarise(hans_m = giw$mean_elevation_m - sum(ele_dis)/sum(inv_dis))
    
  #Create output
  output<-tibble(
    WetID = giw$WetID,
    hans_m = hans$hans_m
  )
  
  #Return output
  output

}

#11.3 Execute function to estimate HANS-----------------------------------------
#Create wrapper function w/ tryCatch for error handling
outside_fun<-function(x){
  tryCatch(fun(x, giws=giws), 
           error = function(e) tibble(WetID = giws$WetID[n],
                                      hans_m = NA))
}

#apply function (~ minutes on SESYNC server)
t0<-Sys.time()
output<-mclapply(X=seq(1,nrow(giws)), FUN=outside_fun, mc.cores=detectCores())
output<-bind_rows(output)
tf<-Sys.time()
tf-t0

#Join to GIWs
giws<-left_join(giws, output)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#12.0 Estimate Network Stats for each wetland===================================
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#12.1 Create network ojbect-----------------------------------------------------
#Create edgelist
edgelist<-giws %>% 
  st_drop_geometry() %>% 
  select(WetID,flow_to) 

#Create network object from edgelist
network<-edgelist %>% graph_from_data_frame()

#12.2 Create function to estimate netowrk order---------------------------------
#Create function to estimate network storage
fun<-function(n){
  
  #isolate wetland
  giw<-giws[n,]
  
  #define upstream paths 
  paths<-all_simple_paths(network, from = paste(giw$WetID), mode = "in")
  paths<-sapply(paths, names) 
  
  #Estiamte longest length
  wet_order<-sapply(paths, length)
  if(length(wet_order)!=0){
    wet_order<-base::max(wet_order, na.rm=T)
    wet_order<-ifelse(wet_order==1, 2, wet_order)
  }else{
    wet_order<-1
  }
  
  #Create simple output
  output<-tibble(
    WetID = giw$WetID, 
    wet_order = wet_order
  )
  
  #Export output
  output
}

#12.3 Estimate wetland order---------------------------------------------------- 
#apply function
output<-mclapply(X=seq(1, nrow(giws)), 
                 FUN = fun, 
                 mc.cores=detectCores()) %>% bind_rows()

#join to giws tibble
giws<-left_join(giws, output)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#13.0 Export Results============================================================
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Export SF shapes
st_write(giws, paste0(data_dir, "III_Products/giws.shp"), delete_layer = T)
st_write(flowpath, paste0(data_dir, "III_Products/flowpath.shp"), delete_layer = T)
st_write(subsheds, paste0(data_dir, "III_Products/subsheds.shp"), delete_layer = T)
st_write(watersheds, paste0(data_dir, "III_Products/watersheds.shp"), delete_layer = T)

#Export raster
writeRaster(dem_burn, paste0(data_dir, "III_Products/dem_burn.tif"), overwrite=T)
save.image('backup_20191103.RData')
