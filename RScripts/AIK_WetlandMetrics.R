#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Title: Wetland Metrics 
#Coder: C. Nathan Jones (cnjones7@ua.edu)
#Date: 6/12/2019
#Purpose: Develop a suite of spatial metrics for AIK's masters thesis
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#To-do list:
#1) Check with Anna to make sure Burn areas make sense
  #Tree island bay    

#1.0 Setup workspace------------------------------------------------------------
#Clear memory
rm(list=ls(all=TRUE))

#Defin relevant working directories
data_dir<-"//nfs/palmer-group-data/choptank/Nate/Spatial_Analysis/"

#Download packages 
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
dem<-raster(paste0(data_dir,"I_Data/2007_1m_DEM.tif"))
  dem[dem<0]<-NA
  dem<-projectRaster(dem, crs=p)  
burn<-st_read(paste0(data_dir,"I_Data/wetland_burn.shp"))
  burn<-st_transform(burn, crs=p)
pnts<-read_csv(paste0(data_dir, "I_Data/q2_locations.csv"))
  pnts<-st_as_sf(pnts, coords  = c('X','Y'), crs = 4269) 
  pnts<-st_transform(pnts, crs=p)
tnc<-st_read(paste0(data_dir,"I_Data/JL_boundary.shp"))
  tnc<-rbind(tnc, st_read(paste0(data_dir,"I_Data/JR_boundary.shp")))
  tnc<-st_transform(tnc, crs=p)
streams<-st_read(paste0(data_dir,"I_Data/TotalStreamNetwork_UCRW.shp"))
  streams<-st_transform(streams, crs=p)
  
#Add Unique ID to points
pnts$WetID<-seq(1,nrow(pnts)) 

#2.0 Identify depressions on TNC Properties-------------------------------------
#Steps
  # (1) Prep Raster
  # (2) Burn AIK sites into DEM
  # (3) Delineate Depressions
  # (4) Filter Depressions

#2.1 Prep raster----------------------------------------------------------------
#Export DEM to workspace
writeRaster(dem, paste0(data_dir,"II_Work/dem.tif"), overwrite=T)

#Fill single cell pits
fill_single_cell_pits(dem=paste0(data_dir,"II_Work/dem.tif"), 
                      output = paste0(data_dir,"II_Work/dem_fill.tif"))

#Conduct Edge preserving filter
edge_preserving_mean_filter(input = paste0(data_dir,"II_Work/dem_fill.tif"), 
                            output = paste0(data_dir,"II_Work/dem_filter.tif"), 
                            threshold = 1)

#2.2 Burn AIK Sites into DEM----------------------------------------------------
#read dem_filter into R memory
dem_filter<-raster(paste0(data_dir,"II_Work/dem_filter.tif"))

#Create burn raster
burn_grd<-fasterize(burn, dem_filter)
burn_grd[burn_grd==1]<-0.5
burn_grd[is.na(burn_grd)]<-1

#Burn into dem
dem_burn<-dem_filter*burn_grd
crs(dem_burn)<-p

#Export burned dem to workspace
writeRaster(dem_burn, paste0(data_dir,"II_Work/dem_burn.tif"), overwrite=T)

#2.3 Delineate "depressions"----------------------------------------------------
#Conduct stochastic depression analysis (This will take ~30 minutes)
stochastic_depression_analysis(dem = paste0(data_dir,"II_Work/dem_burn.tif"), 
                               output = paste0(data_dir,"II_Work/giws.tif"), 
                               rmse = 0.18, 
                               range = 10, 
                               iterations = 100)

#Reclass raster (any depression that was delineated less than 80% of the time is out!)
reclass(input = paste0(data_dir,"II_Work/giws.tif"), 
        output = paste0(data_dir,"II_Work/reclass.tif"), 
        reclass_vals = "'0;0;0.8'")
reclass(input = paste0(data_dir,"II_Work/reclass.tif"), 
        output = paste0(data_dir,"II_Work/reclass.tif"), 
        reclass_vals = "'1;0.8;1'")
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
st_write(giws, paste0(data_dir, "II_Work/giws.shp"), delete_layer=TRUE)

#2.4 Filter delineated depresions-----------------------------------------------
#Filter by area and P:A Ratio
polygon_perimeter(input=paste0(data_dir, "II_Work/giws.shp"))
polygon_area(input=paste0(data_dir, "II_Work/giws.shp"))
giws<-st_read(paste0(data_dir, "II_Work/giws.shp"))
giws<-giws %>%
  #Remove small depressions
  filter(AREA>500) %>%
  #Remove ditched depressions
  mutate(p_a_ratio = AREA/PERIMETER) %>%
  filter(p_a_ratio>3)

#Select depressions within 1 km of TNC properties
giws<-giws[st_buffer(tnc, 1000),]

#Export updated shapefile
st_write(giws, paste0(data_dir, "II_Work/giws.shp"), delete_layer=TRUE)

#3.0 Prepare files for delineation----------------------------------------------
#Steps: 
  # (1) Aggregate to lower resolution DEM
  # (2) Breach depressions
  # (3) Flow direction analysis
  # (4) Flow accumulation analysis

#Aggregate to 10m DEM
aggregate_raster(input = paste0(data_dir,"II_Work/dem_burn.tif"), 
                 output = paste0(data_dir,"II_Work/dem_low_res.tif"),
                 agg_factor = 10)

#Breach Depressions
breach_depressions(dem = paste0(data_dir,"II_Work/dem_low_res.tif"), 
                   output = paste0(data_dir,"II_Work/dem_breach.tif"),
                   fill_pits = T)

#FDR Raster
d8_pointer(dem = paste0(data_dir,"II_Work/dem_breach.tif"), 
           output = paste0(data_dir,"II_Work/fdr_lr.tif"))

#FAC Raster
d8_flow_accumulation(dem = paste0(data_dir,"II_Work/dem_breach.tif"), 
                     output = paste0(data_dir,"II_Work/fac_lr.tif"))

#import into R workspace
fdr<-raster(paste0(data_dir,"II_Work/fdr_lr.tif"))
  crs(fdr)<-p
fac<-raster(paste0(data_dir,"II_Work/fac_lr.tif"))
  crs(fac)<-p
  
#4.0 Delineate Watersheds-------------------------------------------------------
#4.1 Prep GIWs------------------------------------------------------------------
#Add Unique ID for tracking purposes
giws$WetID<-seq(1, nrow(giws))

#4.2 Create function to delineate watershed-------------------------------------
fun<-function(n, #wetland of interest in pnts dataframe
              giws , #wetlands shp derived from DEM [and informed by manual burn shps]){}
              fac , #flow accumulation raster (low res)
              fdr , #flow direction raster (low res)
              dem_hr=dem){ #hi resolution dem
  
  #Steps
  # (1) Create temp directory
  # (2) Delineate watershed with low res DEM
  # (3) Crop high res DEM to "low res" watershed delineation
  # (4) Re-delineate watershed with high res data

  #3.1.1 Create temporary dir~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  scratch_dir<-paste0(tempfile(),"/")
  dir.create(scratch_dir)
  
  #3.1.2 Delineate watershed with low res DEM~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #Identify GIW of interest
  giw<-giws[n,]
  
  #Add 3 m buffer (10 m comes from the resolutin of the raster)
  giw<-st_buffer(giw, 10/3)
  
  #Find max fac point within giw
  fac_giw<-crop(fac, as_Spatial(giw))
  fac_giw<-mask(fac_giw, as_Spatial(giw))
  
  #create pour point
  pp<-rasterToPoints(fac_giw) %>% 
    #convert to tibble
    as_tibble(.) %>%
    #filter to just max fac value
    filter(fac_lr == base::max(fac_lr, na.rm=T)) %>%
    #Select first row
    slice(1)
  
  #Make pour point an sf shape
  pp<-st_as_sf(pp, coords = c("x", "y"), crs = st_crs(giw))
  
  #Export files to scratch dir
  writeRaster(fac, paste0(scratch_dir, "fac.tif"), overwrite=T)
  writeRaster(fdr, paste0(scratch_dir, "fdr.tif"), overwrite=T)
  write_sf(pp, paste0(scratch_dir,"pp.shp"), delete_layer = T)
  
  #Delineate Watershed w/ WBT
  watershed(d8_pntr = paste0(scratch_dir, "fdr.tif"), 
            pour_pts = paste0(scratch_dir,"pp.shp"), 
            output = paste0(scratch_dir,"watershed.tif"))
  
  #Read watershed in
  w_grd<-raster(paste0(scratch_dir,"watershed.tif"))
  
  #Convert watershed to polygon
  w_shp<- w_grd %>% st_as_stars() %>% st_as_sf(., merge = TRUE)
  
  #3.1.3 Re-delineate watershed with high res data~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #Clip high resolution DEM with w_shp
  mask<-st_buffer(w_shp, 100)
  dem_hr<-crop(dem_hr, as_Spatial(mask))
  
  #Write raster to scratch workspace
  writeRaster(dem_hr,paste0(scratch_dir,"dem_hr.tif"), overwrite=T)
  
  #Breach Analysis of the DEM
  breach_depressions(dem = paste0(scratch_dir,"dem_hr.tif"), 
                     output = paste0(scratch_dir,"dem_breach.tif"), 
                     fill_pits = TRUE)
  
  #Flow Direction and flow accumulation
  d8_pointer(dem = paste0(scratch_dir,"dem_breach.tif"), 
             output = paste0(scratch_dir,"fdr_hr.tif"))
  d8_flow_accumulation(dem = paste0(scratch_dir,"dem_breach.tif"), 
                       output = paste0(scratch_dir,"fac_hr.tif"))
  
  #Read fac and fdr into R environment
  fac_hr<-raster(paste0(scratch_dir,"fac_hr.tif"))
  fdr_hr<-raster(paste0(scratch_dir,"fdr_hr.tif"))
  
  #Find max fac point within giw
  fac_giw_hr<-crop(fac_hr, as_Spatial(giw))
  fac_giw_hr<-mask(fac_giw_hr, as_Spatial(giw))
  
  #create pour point
  pp_hr<-rasterToPoints(fac_giw_hr) %>% 
    #convert to tibble
    as_tibble() %>%
    #select point with max fac
    filter(fac_hr==base::max(fac_hr))
  
  #Make pour point an sf shape
  pp_hr<-st_as_sf(pp_hr, coords = c("x", "y"), crs = st_crs(giw))
  
  #Export pour point to scratch workspace
  write_sf(pp_hr, paste0(scratch_dir,"pp_hr.shp"), delete_layer = T)
  
  #Delineate Watershed w/ WBT
  watershed(d8_pntr = paste0(scratch_dir, "fdr_hr.tif"), 
            pour_pts = paste0(scratch_dir,"pp_hr.shp"), 
            output = paste0(scratch_dir,"watershed_hr.tif"))
  
  #Read watershed in
  w_hr_grd<-raster(paste0(scratch_dir,"watershed_hr.tif"))
  
  #Convert watershed to polygon
  w_hr_shp<- w_hr_grd %>% st_as_stars() %>% st_as_sf(., merge = TRUE)
  
  #If there are multiple shapes, then combine
  if(nrow(w_hr_shp)>1){
    w_hr_shp<-st_union(w_hr_shp)
  }
  
  #Add id info and write to output folder
  write_sf(w_hr_shp, paste0(data_dir, "II_Work/",giw$WetID,"_high_res_watershed.shp"), delete_layer = T)
  
  #Export watershed size
  data.frame(WetID=giw$WetID, ws_area=paste(st_area(w_hr_shp)))
}

#4.3 Apply function-------------------------------------------------------------
#Create wrapper function w/ tryCatch for error handling
outside_fun<-function(x){
  tryCatch(fun(x, giws=giws, fac=fac,fdr=fdr,dem_hr=dem), 
           error = function(e) data.frame(WetID=x, ws_area="-999"))
}
  
#apply function
output<-mclapply(X=seq(1,nrow(giws)), FUN=outside_fun, mc.cores=detectCores())
output<-bind_rows(output)

#Join to master GIW tibble
giws<-left_join(giws, output)

#Gather shapes and export~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Make list of watershed shape files
watersheds<-list.files(paste0(data_dir,"II_work")) %>%
  enframe() %>%
  select(value) %>%
  filter(str_detect(value, "high_res_watershed.shp"))

#load those shapes
watersheds_shp<-st_read(paste0(data_dir,"II_work/", watersheds[1,]))
watersheds_shp$ID<-watersheds$value[1] %>% str_match_all("[0-9]+") %>% unlist
jig<-colnames(watersheds_shp)
for(i in 2:nrow(watersheds)){
  print(i)
  temp_shp<-st_read(paste0(data_dir,"II_work/", watersheds[i,]))
  temp_shp$ID<-watersheds$value[i] %>% str_match_all("[0-9]+") %>% unlist
  colnames(temp_shp)<-jig
  watersheds_shp<-rbind(watersheds_shp, temp_shp)
}

#export 
st_write(watersheds_shp, paste0(data_dir,"III_Products/watersheds.shp"), delete_layer =T)

#5.0 Wetland Subshed Delineation------------------------------------------------
#5.1 Create function to identify indididual subsheds----------------------------
fun<-function(n,dem){
#Steps
  # (1) Identify watershed
  # (2) Select wetlands within the watershed
  # (3) Load watersheds draining to those wetlands
  # (4) Substract using raster algebra magic

  #Identify giw
  giw<-giws[n,]
  
  #Download watershed shape
  ws_shp<-st_read(paste0(data_dir, "II_Work/",giw$WetID, "_high_res_watershed.shp"))
  ws_grd<-fasterize(ws_shp, dem)
  ws_grd<-crop(ws_grd, ws_shp)
  
  #Identify upslope depressions
  giw_upslope<-giws[ws_shp,] %>% filter(WetID!=giw$WetID)
  
  #if there are upslope wetlands, then substract their subsheds!!!!
  if(nrow(giw_upslope)>0){
    #Iteratively substract uplsope subshed areas
    for(i in 1:nrow(giw_upslope)){
      #Identify upslope subshed
      temp_shp<-st_read(paste0(data_dir, "II_Work/",giw_upslope$WetID[i], "_high_res_watershed.shp"))
      temp_grd<-fasterize(temp_shp, ws_grd)
      
      #Substract upslope subshed from larger watershed
      temp_grd[temp_grd==1]<-0
      temp_grd[is.na(temp_grd)]<-1
      ws_grd<-ws_grd*temp_grd
    }
    
    #Convert subshed grd to shp
    ws_grd[ws_grd==0]<-NA
    ws_shp<-ws_grd %>% st_as_stars() %>% st_as_sf(., merge = TRUE)
    
    #Just incase there are multiple shapes
    ws_shp<-ws_shp[giw,]
  }
  
  #Add id info and write to output folder
  write_sf(ws_shp, paste0(data_dir, "II_Work/",giw$WetID,"_subshed.shp"), delete_layer = T)
  
  #Export watershed size
  data.frame(WetID=giw$WetID, subshed_area=paste(st_area(ws_shp)))
}

#5.2 Apply function-------------------------------------------------------------
#Create wrapper function w/ tryCatch for error handling
outside_fun<-function(x){
  tryCatch(fun(x, dem=dem), 
           error = function(e) data.frame(WetID=x, subshed_area="-9999"))
}

#apply function (~2 minutes on SESYNC server)
output<-mclapply(X=seq(1,nrow(giws)), FUN=outside_fun, mc.cores=detectCores())
output<-bind_rows(output) %>%
  group_by(WetID) %>%
  summarise(subshed_area = base::max(subshed_area, na.rm=T))

#Merge with GIW tible
giws<-left_join(giws, output)

#5.3 Gather Shapes--------------------------------------------------------------
#Gather shapes and export~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Make list of watershed shape files
subsheds<-list.files(paste0(data_dir,"II_work")) %>%
  enframe() %>%
  select(value) %>%
  filter(str_detect(value, "subshed.shp"))

#load those shapes
subsheds_shp<-st_read(paste0(data_dir,"II_work/", subsheds[1,]))
subsheds_shp$ID<-subsheds$value[1] %>% str_match_all("[0-9]+") %>% unlist
colnames(subsheds_shp)<-c("shape","geometry","ID")
for(i in 2:nrow(subsheds)){
  print(i)
  temp_shp<-st_read(paste0(data_dir,"II_work/", subsheds[i,]))
  if(nrow(temp_shp)>0){
    temp_shp$ID<-subsheds$value[i] %>% str_match_all("[0-9]+") %>% unlist
    colnames(temp_shp)<-c("shape","geometry","ID")
    subsheds_shp<-rbind(subsheds_shp, temp_shp)
  }
}

#export 
st_write(subsheds_shp, paste0(data_dir,"III_Products/subsheds.shp"), delete_layer =T)


#6.0 Calculate Storage Capacity-------------------------------------------------
#6.1 Create function to estimate storage capacity ------------------------------
fun<-function(n, dem, dz, z_max){
  #Steps
  # (1) Gather GIW, DEM, and subshed data
  # (2) Crop DEM + Locate minum elevation in depression
  # (4) Inundation function
  # (5) Execute and find break point
  # (6) Print inundation
  # (7) Export!
  
  #Identify wetland of interest
  giw<-giws[n,]

  #Call subshed shape
  w_shp<-st_read(paste0(data_dir, "II_Work/",giw$WetID,"_subshed.shp"))

  #Crop DEM to subshed
  temp<-crop(dem, w_shp)
  temp<-mask(temp, w_shp)
  
  #Create Minimum Raster
  temp_min<-temp*0+minValue(temp)
  temp_min@crs<-dem@crs
  
  #Mark deepest point in delineated wetland
  pnt_min<-crop(temp, giw)
  pnt_min<-mask(pnt_min, giw)
  pnt_min<-rasterToPoints(pnt_min, fun = function(x) x==cellStats(pnt_min, base::min)) 
  pnt_min<-data.frame(x=pnt_min[1,1],y=pnt_min[1,2])
  pnt_min<-st_as_sf(pnt_min, coords = c("x", "y"), crs = st_crs(giw))
  
  #Create function to return conditional raster
  Con<-function(condition, trueValue, falseValue){
    return(condition * trueValue + (!condition)*falseValue)
  }
  
  #Create function to calcluate inundation area, volume, and spill boundary length
  inundate<-function(z){
    
    #define inundated area [area connected to pnt_min]
    area<-Con(temp>(temp_min+z),0,1)
    area<-raster::clump(area)
    group<-raster::extract(x=area, y=pnt_min)
    area[area!=group]<-NA
    area[area==group]<-1
    
    #Create metrics to estimate area, volume, and spill boundary length
    volume<-(((z+temp)-temp_min)*area)*res(area)[1]*res(area)[2]
    outflow<-cellStats(area*boundaries(temp_min, type="inner"), 'sum')
    
    #Export Data
    c(z, #inundation depth
      cellStats(area, 'sum')*res(area)[1]*res(area)[2], #area (m^2)
      cellStats(volume, 'sum'), #volume (m^3)
      outflow) #Outflow length (3 m increments)
  }
  
  #Create function to calculate inundation area/volume
  df<-lapply(seq(dz,z_max,dz),inundate)
  df<-do.call(rbind, df)
  df<-data.frame(df)  
  colnames(df)<-c("z", "area","volume","outflow_length")
  
  #Identify spillpoint
  df<-df %>% filter(outflow_length>0) %>% slice(1)
  
  #Print max inundation
  area<-Con(temp>(temp_min+df$z),0,1)
  area<-raster::clump(area)
  group<-raster::extract(x=area, y=pnt_min)
  area[area!=group]<-NA
  area[area==group]<-1
  area<-area %>% st_as_stars() %>% st_as_sf(., merge = TRUE)
  st_write(area, paste0(data_dir, "II_Work/",giw$WetID,"_inundation.shp"), delete_layer = T)
  
  #Assign ID
  df<- df %>%
    mutate(WetID = giw$WetID) %>%
    mutate(min_ele = cellStats(temp_min, base::min)) %>%
    rename(inun_area = area)
  
  #Export inundation information
  df
}
  
#6.2 Apply function-------------------------------------------------------------
#Create wrapper function w/ tryCatch for error handling
outside_fun<-function(x){
  tryCatch(fun(x, dem=dem, dz = 0.05, z_max=2), 
           error = function(e) data.frame(z = -9999, 
                                          inun_area = -9999, 
                                          volume = -9999,
                                          outflow_length = -9999, 
                                          WetID=x))
}

#apply function (~ minutes on SESYNC server)
t0<-Sys.time()
output<-mclapply(X=seq(1,nrow(giws)), FUN=outside_fun, mc.cores=detectCores())
output<-bind_rows(output)
tf<-Sys.time()
tf-t0

#Join to GIWs
colnames(output)<-c("z_max", "inun_area", "inun_volume", "outflow_length", "WetID", "min_ele")
giws<-left_join(giws, output)

#Save image as backup 
save.image("backup.RData")

#7.0 Estimate the distribution of storage capacity------------------------------
#7.1 Create function------------------------------------------------------------
fun<-function(
  n,
  giws, 
  watersheds_shp){

  #Filter -9999 values from giw tibble
  giws<-giws %>% filter(inun_volume>0)
  
  #Gather relevant data
  giw<-giws[n,]
  subshed<-st_read(paste0(data_dir, "II_Work/",giw$WetID,"_subshed.shp"))
  watershed<-st_read(paste0(data_dir, "II_Work/",giw$WetID,"_high_res_watershed.shp"))
  
  #Identify giws in watershed
  siblings<-giws[watershed,] %>% filter(WetID!=giw$WetID)
  
  #Identify neighborhood [for now, we're defining the neighborhood as adjacent watersheds]
  neighborhood<-watersheds_shp[st_buffer(watershed, 10),]
  neighbors<-giws[neighborhood,] %>% filter(WetID!=giw$WetID)
  
  #Estiamte storage capacity [units = cm] 
  s_subshed<-giw$inun_volume/as.numeric(giw$subshed_area)*100
  s_watershed<-sum(siblings$inun_volume, na.rm=T)/as.numeric(st_area(watershed))*100
  s_neighborhood<-sum(neighbors$inun_volume, na.rm=T)/sum(as.numeric(st_area(neighborhood)), na.rm=T)*100
  
  #Estimate volume ratio
  s_ratio_ws<-giw$inun_volume/sum(siblings$inun_volume, na.rm=T)
  s_ratio_nh<-sum(siblings$inun_volume, na.rm=T)/sum(neighbors$inun_volume, na.rm=T)
  
  #Estimate vertical difference (units = cm)
  s_vert_ws<-(mean(siblings$min_ele, na.rm=T) - giw$min_ele)*100
  s_vert_nh<-(mean(neighbors$min_ele, na.rm=T) - giw$min_ele)*100
  
  #gather results
  t<-tibble(
    WetID = giw$WetID,
    s_subshed, 
    s_watershed,
    s_neighborhood, 
    s_ratio_ws, 
    s_ratio_nh,
    s_vert_ws, 
    s_vert_nh
  )
  
  #Export
  t
}

#7.2 Execute function-----------------------------------------------------------
#Create wrapper function w/ tryCatch for error handling
outside_fun<-function(x){
  tryCatch(fun(x, giws=giws, watersheds_shp = watersheds_shp), 
           error = function(e) tibble(WetID = giws$WetID[n],
                                      s_subshed = NA, 
                                      s_watershed = NA,
                                      s_neighborhood = NA, 
                                      s_ratio_ws = NA, 
                                      s_ratio_nh = NA,
                                      s_vert_ws = NA, 
                                      s_vert_nh = NA))
}

#apply function (~ minutes on SESYNC server)
t0<-Sys.time()
output<-mclapply(X=seq(1,nrow(giws)), FUN=outside_fun, mc.cores=detectCores())
output<-bind_rows(output)
tf<-Sys.time()
tf-t0

#Join to GIWs
giws<-left_join(giws, output)
