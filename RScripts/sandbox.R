#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Title: Sandbox
#Coder: C. Nathan Jones (cnjones7@ua.edu)
#Date: 8/14/2019
#Purpose: Examine 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#To-do
#1) make subshed delineation based on unbreached dem [maybe?]

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
streams<-st_read(paste0(data_dir, "I_Data/streams_jr.shp"))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#2.0 Filter the DEM=============================================================
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#2.1 Fill Single Cell Pitts-----------------------------------------------------
#Export DEM to workspace
writeRaster(dem, paste0(data_dir,"II_Work/dem.tif"), overwrite=T)

#Fill single cell pits
fill_single_cell_pits(dem=paste0(data_dir,"II_Work/dem.tif"), 
                      output = paste0(data_dir,"II_Work/dem_fill.tif"))

#Read raster back into workspace
dem_fill<-raster(paste0(data_dir,"II_Work/dem_fill.tif"))

#2.2 Filter DEM-----------------------------------------------------------------
#Apply simple gausian filter to smooth random errors from DEM
dem_filter<- focal(dem_fill, w=focalWeight(dem, 3, "Gauss"))
crs(dem_filter)<-p

#2.3 Burn Streams into DEM------------------------------------------------------
#Convert streams to single polyline
streams_grd<-fasterize(streams, dem_filter)

#Create upland mask
upland_grd<-streams_grd*0
upland_grd[is.na(upland_grd)]<-1
upland_grd[upland_grd==0]<-NA

#Add addition raster to DEM
dem_upland<- dem_filter*upland_grd

#Fill dem_add
writeRaster(dem_upland, paste0(data_dir, "II_Work/dem_upland.tif"), overwrite=T)
fill_missing_data(input = paste0(data_dir, "II_Work/dem_upland.tif"), 
                  output = paste0(data_dir, "II_Work/dem_burn.tif"))
dem_burn<-raster(paste0(data_dir, "II_Work/dem_burn.tif"))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#3.0 Define "Root" Depressions==================================================
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#3.1 Use the Stochastic Deprresion Tool to identify deprresions-----------------
#Export fitlered DEM to workspace
#writeRaster(dem_burn, paste0(data_dir,"II_Work/dem_burn.tif"), overwrite=T)

#Apply stochastic depressin analysis tool
set.seed(100)
stochastic_depression_analysis(dem = paste0(data_dir,"II_Work/dem_burn.tif"), 
                               output = paste0(data_dir,"II_Work/giws.tif"), 
                               rmse = 0.18, 
                               range = 10, 
                               iterations = 1000)

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
fun<-function(n, giws, dem, max_size = 250, n_slices = 10){
  
  #Identify GIW
  giw<-giws[n,]
  
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
  tryCatch(fun(x, giws, dem, max_size = 250, n_slices = 10), 
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
writeRaster(dem_burn, paste0(data_dir, "dem.tif"), overwrite=T)

#Execute breach tool
breach_depressions(dem = paste0(data_dir,"dem.tif"), 
                   output = paste0(data_dir,"dem_breach.tif"), 
                   fill_pits = TRUE)

#Execute WBT flow accumulation raster
d8_flow_accumulation(dem = paste0(data_dir,"dem_breach.tif"), 
                     output = paste0(data_dir,"fac.tif"))

d8_pointer(dem = paste0(data_dir, "dem_breach.tif"), 
           output = paste0(data_dir, "fdr.tif"))

#Pull fac back into R environment
fac<-raster(paste0(data_dir, "fac.tif"))
fdr<-raster(paste0(data_dir, "fdr.tif"))

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
    giw_bndry<- giws_ds %>% select(WetID) %>% st_cast(., "POINT")
    
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
st_write(seeds, paste0(data_dir, "seeds.shp"), delete_layer = TRUE)

#Execute WBT flowpaths function
trace_downslope_flowpaths(seed_pts = paste0(data_dir,"seeds.shp"), 
                          d8_pntr  = paste0(data_dir, "fdr.tif"), 
                          output   = paste0(data_dir, "flowpath.tif"))

#Read flowpath raster into R environment
flowpath<-raster(paste0(data_dir, "flowpath.tif"))

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
st_write(seeds, paste0(data_dir, "seeds.shp"), delete_layer = TRUE)

#Execute WBT flowpaths function
trace_downslope_flowpaths(seed_pts = paste0(data_dir,"seeds.shp"), 
                          d8_pntr  = paste0(data_dir, "fdr.tif"), 
                          output   = paste0(data_dir, "flowpath.tif"))

#Read flowpath raster into R environment
fp<-raster(paste0(data_dir, "flowpath.tif"))*0+1

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
st_write(pp, paste0(data_dir, "pp.shp"), delete_layer=T)
writeRaster(dem_burn, paste0(data_dir, "dem.tif"), overwrite=T)

#Conduct watershed analysis
watershed(d8_pntr = paste0(data_dir, "fdr.tif"), 
          pour_pts = paste0(data_dir,"pp.shp"), 
          output = paste0(data_dir,"subshed.tif"))

#Read subsheds into R environment
subsheds<-raster(paste0(data_dir,"subshed.tif"))

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
  
  #Use igraph to determine upstream subsheds
  edgelist<-giws %>% 
    st_drop_geometry() %>% 
    select(WetID,flow_to) %>% 
    mutate(flow_to = if_else(is.na(flow_to), 
                             0,
                             flow_to))
  
  network<-edgelist %>% graph_from_data_frame()
  paths_in<-all_simple_paths(network, from = UID, mode = "in")
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
output<-mclapply(X=seq(1,nrow(giws)), FUN=outside_fun, mc.cores=detectCores()) %>% do.call(rbind,.)


#Add WetID to subsheds
watersheds<-left_join(watersheds, output)


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

#7.3.3 Watersheds~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Use a spatial join function for now [graduate to network analysis later?]

#Creat function to find wetlands within a given wateshed
#fun<-function(n, giws, watersheds)

#Select of interest 
giw<-giws %>% slice(n)

#Define its WetID [for later use]
WetID<-giw$WetID

#Identify root wetland
#Is this a root wetlands?
if(giw$WetID != giw$root_giw){
  #If not, redefine giw as its root wetland
  giw<-giws %>% filter(WetID == giw$root_giw)
}

#Identify watershed
watershed<-watersheds %>% filter(WetID == giw$WetID)

#Welp...that's a good place to quit.  My watershed IDs and WetIDs are a mis-match. Soo..
#  tomorrow: (1) Fix that shit, 
#            (2) Identify all wetlands withint watershed by their centroids
#            (3) Sum everything up adn divide by the ws_area!
#            (4) Remember, when you spit it out, use WetID defined above as ID...
#  then, calculate vertical metrics and send to Anna.

#Identify matching watershed

#Create centroid of giws
centroids<-giws %>% st_centroid(., by_element=T) %>% select(WetID)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Plot for funzies===============================================================
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
plot(dem)
streams %>% plot(., col="dodgerblue4", add=T)
giws %>% filter(is.na(merge_to)) %>% st_geometry() %>% plot(., add=T, lty=2, col="dodgerblue4")
giws %>% filter(spawned == 0)  %>% st_geometry() %>% plot(., add=T, col="dodgerblue2")
subsheds %>% st_geometry() %>% plot(., add=T, lcol="grey80", lwd=0.5)
flowpath %>% st_geometry() %>% plot(., add=T, border="dodgerblue2")
pp %>% st_geometry() %>% plot(., add=T, pch=19, col="grey30")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Code Graveyard=================================================================
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#6.3.1 Watershed Delineation~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Create a temporary workspace
scratch_dir<-paste0(tempfile(),"/")
dir.create(scratch_dir)

#Run WBT unnested watershed delineation tool
unnest_basins(d8_pntr = paste0(data_dir, "fdr.tif"), 
              pour_pts = paste0(data_dir,"pp.shp"), 
              output = paste0(scratch_dir,"watershed.tif"))

#Create list of raster files
files<-list.files(scratch_dir)

#Read basins and convert to polygons
watersheds<-raster(paste0(scratch_dir,files[1])) %>% 
  st_as_stars(.) %>% st_as_sf(., merge=T) %>%
  rename(ws_id = substr(files[1],1, nchar(files[1])-4))
for(i in 2:length(files)){
  print(i)
  temp<-raster(paste0(scratch_dir,files[i])) %>% 
    st_as_stars(.) %>% st_as_sf(., merge=T) %>%
    rename(ws_id = substr(files[i],1, nchar(files[i])-4))
  watersheds<-rbind(watersheds, temp)
}

#Delete temporary workspace
unlink(scratch_dir, recursive = T)

#6.3.2 Add WetID to Watersheds~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Create function to identify subhsed
fun<-function(n, 
              watersheds,
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
  
  #Identify overlapping watershed
  watershed<-watersheds[pp_new,]
  
  #Create output of subhsed id and WetID
  output<-tibble(
    WetID = p$WetID, 
    ws_id = watershed$ws_id
  )
  
  #Export output
  output
}

#Apply function
outside_fun<-function(x){
  tryCatch(fun(x, watersheds, pp, fp), 
           error = function(e) NA)
}

#apply function (~ minutes on SESYNC server)
t0<-Sys.time()
output<-mclapply(X=seq(1,nrow(pp)), FUN=outside_fun, mc.cores=detectCores())
output<-bind_rows(output)
tf<-Sys.time()
tf-t0

#Add WetID to subsheds
watersheds<-left_join(watersheds, output)

