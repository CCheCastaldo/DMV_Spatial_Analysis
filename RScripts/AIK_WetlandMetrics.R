#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Title: Wetland Metrics 
#Coder: C. Nathan Jones (cnjones7@ua.edu)
#Date: 6/12/2019
#Purpose: Develop a suite of spatial metrics for AIK's masters thesis
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#To-do list:
#0) Update to latest WBT exe
#1) Use 2007 DEM to be consistent with KH
#2) Check with Anna to make sure Burn areas make sense
#3) Convert from sp hell to sf
#4) Convert to metric units

#1.0 Setup workspace------------------------------------------------------------
#Clear memory
rm(list=ls(all=TRUE))

#Defin relevant working directories
data_dir<-"//nfs/palmer-group-data/choptank/Nate/Storage_Capacity/"

#Download packages 
library(tidyverse)
library(raster)
library(sf)
library(parallel)
library(whitebox)
library(fasterize)
library(stars)

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
burn<-st_read(paste0(data_dir,"II_Work/wetland_burn.shp"))
  burn<-st_transform(burn, crs=p)
pnts<-read_csv(paste0(data_dir, "I_Data/q2_locations.csv"))
  pnts<-st_as_sf(pnts, coords  = c('X','Y'), crs = 4269) 
  pnts<-st_transform(pnts, crs=p)

#Add Unique ID to points
pnts$WetID<-seq(1,nrow(pnts)) 

#2.0 Prepare files for delineation----------------------------------------------
#Steps: 
  #(1) Filter DEM
  #(2) Burn Ponds into Raster
  #(3) Breach depressions
  #(4) Flow direction analysis
  #(5) Flow accumulation analysis

#2.1 Create low resolution dem--------------------------------------------------
writeRaster(dem, paste0(data_dir,"II_Work/dem.tif"), overwrite=T)
aggregate_raster(input = paste0(data_dir,"II_Work/dem.tif"), 
                 output = paste0(data_dir,"II_Work/dem_10m.tif"), 
                 agg_factor = 10)
dem_lr<-raster(paste0(data_dir,"II_Work/dem_10m.tif"))
dem_lr[dem_lr<0]<-NA

#2.2 Burn known wetlands into DEM-----------------------------------------------
#Create burn raster
burn_grd<-fasterize(burn, dem_lr)
  burn_grd[burn_grd==1]<-0
  burn_grd[is.na(burn_grd)]<-1

#Burn into dem
dem_burn<-dem_lr*burn_grd
  crs(dem_burn)<-p

#Write burn to scratch workspace
writeRaster(dem_burn, paste0(data_dir,"II_Work/dem_burn.tif"), overwrite=T)

#2.3 Define depressions (again, low res)----------------------------------------
#Use the stochastic depression analysis tool!
stochastic_depression_analysis(dem = paste0(data_dir,"II_Work/dem_burn.tif"), 
                               output = paste0(data_dir,"II_Work/giws_lr.tif"), 
                               rmse = 0.18, 
                               range = 10, 
                               iterations = 100)

#Reclass raster (any depression that was delineated less than 80% of the time is out!)
reclass(input = paste0(data_dir,"II_Work/giws_lr.tif"), 
        output = paste0(data_dir,"II_Work/reclass.tif"), 
        reclass_vals = "'0;0;0.8'")
reclass(input = paste0(data_dir,"II_Work/reclass.tif"), 
        output = paste0(data_dir,"II_Work/reclass.tif"), 
        reclass_vals = "'1;0.8;1'")
r<-raster(paste0(data_dir,"II_Work/reclass.tif")) 
r[r==0]<-NA
writeRaster(r, paste0(data_dir, "II_Work/reclass_na.tif"), overwrite=T)
remove(r)

#Identify "clumps" of depressions
clump(input = paste0(data_dir, "II_Work/reclass_na.tif"), 
      output = paste0(data_dir, "II_Work/group.tif"), 
      zero_back = T)

#Identify Clusters greater than 100m2
c<-raster(paste0(data_dir, "II_Work/group.tif")) 
c[c==1]<-NA
c<- c %>% 
  #Convert to polygon
  st_as_stars(.) %>% st_as_sf(., merge = TRUE) %>%
  #Estiamte polygon areas
  mutate(area = as.numeric(st_area(., by_id=T))) %>%
  #filter small areas out
  filter(area>100) 
st_write(c, paste0(data_dir, "II_Work/giws_lr.shp"), delete_layer=TRUE)
remove(c)

#2.4 Create low resolution fdr and fac------------------------------------------
#Breach depressions
breach_depressions(dem = paste0(data_dir,"II_Work/dem_burn.tif"), 
                   output = paste0(data_dir,"II_Work/dem_breach.tif"))

#FDR Raster
d8_pointer(dem = paste0(data_dir,"II_Work/dem_breach.tif"), 
           output = paste0(data_dir,"II_Work/fdr_lr.tif"))

#FAC Raster
d8_flow_accumulation(dem = paste0(data_dir,"II_Work/dem_breach.tif"), 
                     output = paste0(data_dir,"II_Work/fac_lr.tif"))

#import into R workspace
burn<-raster(paste0(data_dir,"II_Work/dem_burn.tif"))
breach<-raster(paste0(data_dir,"II_Work/dem_breach.tif"))
fdr<-raster(paste0(data_dir,"II_Work/fdr_lr.tif"))
fac<-raster(paste0(data_dir,"II_Work/fac_lr.tif"))

#3.0 Delineate Watersheds-------------------------------------------------------
#Create function to delineate watershed
#fun<-function(){}

#Create temporary file
temp_dir<-paste0(tempfile(),"/")
dir.create(temp_dir)










#3.0 Delineate Wetlands---------------------------------------------------------
#Create function to identify depressions
fun<-function(dem=dem,
              min_size=5000, #minimum area (in map units)
              iterations=100, #number of Monte Carlo iterations
              dem_rmse=0.0607){ #RMSE of 18.5 cm

#Create temporary file
temp_dir<-paste0(tempfile(),"/")
dir.create(temp_dir)
  
#Write dem to scratch workspace
writeRaster(dem, paste0(temp_dir, "dem.tif"), overwrite=T)

#Fill single-cell depressions
fill_single_cell_pits(dem    = paste0(temp_dir, "dem.tif"), 
                      output = paste0(temp_dir,"dem_fill.tif"))

#Filter the DEM
edge_preserving_mean_filter(input  = paste0(temp_dir,"dem_fill.tif"), 
                            output = paste0(temp_dir,"dem_filter.tif"), 
                            threshold = 100)

#Identify depressions
stochastic_depression_analysis(dem = paste0(temp_dir,"dem_filter.tif"), 
                               output = paste0(temp_dir,"depression.tif"), 
                               rmse = dem_rmse,
                               range = 100, 
                               iterations =iterations)

#Reclass raster (any depression that was delineated less than 80% of the time is out!)
reclass(input = paste0(temp_dir,"depression.tif"), 
        output = paste0(temp_dir,"reclass.tif"), 
        reclass_vals = "'0;0;0.8'")
reclass(input = paste0(temp_dir,"reclass.tif"), 
        output = paste0(temp_dir,"reclass.tif"), 
        reclass_vals = "'1;0.8;1'")
r<-raster(paste0(temp_dir,"reclass.tif")) 
r[r==0]<-NA
writeRaster(r, paste0(temp_dir, "reclass_na.tif"), overwrite=T)

#Identify "clumps" of depressions
clump(input = paste0(temp_dir,"reclass_na.tif"), 
      output = paste0(temp_dir,"group.tif"), 
      zero_back = T)

#Identify Clusters greater than 100m2
w<-raster(paste0(temp_dir,"group.tif")) 
w[w==1]<-NA
w<- w %>% 
  #Convert to polygon
  st_as_stars(.) %>% st_as_sf(., merge = TRUE) %>%
  #Estiamte polygon areas
  mutate(area = as.numeric(st_area(., by_id=T))) %>%
  #filter small areas out
  filter(area>min_size) 

#Export
w
}

#Apply function 
giws_jl<-fun(dem=dem_burn_jl)
giws_jr<-fun(dem=dem_burn_jr)

#combine giws
giws<-rbind(giws_jl, giws_jr)

#export to data directory
st_write(giws_jl, paste0(data_dir, "II_Work/giws_jl.shp"))
st_write(giws_jr, paste0(data_dir, "II_Work/giws_jr.shp"))
st_write(giws,    paste0(data_dir, "II_Work/giws.shp"))

#4.0 Delinate Watersheds--------------------------------------------------------
#Create function to delineate watershed
#fun<-function(n, 
#              pnts = pnts, 
#              giws = giws){}

#For testing
n<-1

#Isolate wetland of interests
giw<-giws[pnts[n,],]


















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

#4.0 Calculate Storage Capacity-------------------------------------------------
#4.1 Stage storage relationships------------------------------------------------
fun<-function(WetID){
  #Define point of interest
  pnt<-pnts[pnts$WetID==WetID,]
  
  #Define dem and depressions of interest
  if(gIntersects(pnt, burn_jr)==T){
    dem=dem_jr
  }else{
    dem=dem_jl
  }
  
  #Define subshed 
  subshed<-raster(paste0(data_dir, "II_work/subshed_",WetID,".asc"))
  
  #Devleop stage storage relationship
  storage<-GIW_stage_storage(
    subshed = subshed, #Subshed delineated in previous step
    dem = dem,         #DEM
    z_max = 2,         #Maximum inundation depth
    dz = 1/12)         #Inundation Interval
  
  #Export
  storage$WetID<-WetID
  storage
}

#Run function
t0<-Sys.time()
cl <- makePSOCKcluster(detectCores()) #Create Clusters
clusterEvalQ(cl, library(raster))   
clusterEvalQ(cl, library(rgeos))
clusterExport(cl, c('pnts','burn_jr','dem_jr',"dem_jl", "data_dir",'GIW_stage_storage'), env=.GlobalEnv)  #Send Clusters function with the execute function
x<-parLapply(cl, seq(1,nrow(burn)), fun) #Run execute Function
stopCluster(cl)  #Turn clusters off
tf<-Sys.time()
tf-t0

#unlist output
storage<-do.call(rbind,x)

#4.2 Spill point----------------------------------------------------------------
#Create function to estimate the spill point and print plot 
fun<-function(WetID){
  #Define point of interest
  pnt<-pnts[pnts$WetID==WetID,]
  
  #Define dem and depressions of interest
  if(gIntersects(pnt, burn_jr)==T){
    dem=dem_burn_jr
  }else{
    dem=dem_burn_jl
  }
  
  #Define subshed 
  subshed<-raster(paste0(data_dir, "II_work/subshed_",WetID,".asc"))
  #Define Storage Curve
  storage<-storage[storage$WetID==WetID,]
  
  #Calculate max inundation
  inundation<-GIW_max_inundation(
    subshed,  #wateshed raster
    dem, #DEM for the analysis
    storage   #Storage Curve 
  )
  
  #Create layers for plotting
  subshed<-rasterToPolygons(subshed, dissolve=T)
  dem<-crop(dem, subshed)
  dem<-mask(dem, subshed)  
  
  #Estimate Below water surface (using W and Lane 2016)
  a_ha<-(storage$area[1]*(.3048^2))*0.0001
  dV_ham=0.25*(a_ha^1.4742)
  dV<-(dV_ham*10000)*(3.28084^3)
  storage$volume<-storage$volume+dV
  
  #Create Output Plots!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  png(paste0(scratch_dir,WetID,".png"), 
      width = 6, height=3, units="in", res=150)
  
  #Plot side by side
  par(mfrow=c(1,2))
  
  #Plot Inundation
  par(mar=c(1,1,1,1))
  plot(subshed, main=paste0("Wetland = ",pnt$wetland))
  plot(dem, add=T, legend=F)
  plot(inundation, add=T, legend=F, col="blue")
  
  #Plot stage storage curve
  par(mar=c(4,4,1,1))
  par(mgp=c(2.1,1,0))
  plot(storage$z*12, storage$volume/gArea(subshed)*12, type="n",
       ps=12, cex.axis=10/12, cex.lab=14/12, 
       xlab="Relative Wetland Stage [in]", 
       ylab=expression("Specific Storage [in"^3*"/in"^2*"]"))
  abline(h=storage$volume[storage$outflow_length>0][1]/gArea(subshed)*12, 
         lty=2, lwd=2, col="grey30")
  points(storage$z*12, storage$volume/gArea(subshed)*12,
         pch=21, col="grey30", bg="grey70", cex=2)
  
  #Turn device off
  dev.off()
  
  #Create df of relevant storage capacity information 
  output<-data.frame(WetID=WetID,
                     area_subshed_ft2=gArea(subshed),
                     area_wetland_ft2=storage$area[storage$outflow_length>0][1], 
                     volume_ft3=storage$volume[storage$outflow_length>0][1])
  #Export!
  output
}


#Run function
output<-lapply(seq(1,nrow(burn)), fun)
output<-do.call(rbind,output)

#5.0 Export results-------------------------------------------------------------
#Combine points with output data
pnts<-pnts@data
pnts<-left_join(pnts, output)

#Export CSV file
write_csv(pnts, paste0(data_dir, "III_products/AIK_output.csv"))
