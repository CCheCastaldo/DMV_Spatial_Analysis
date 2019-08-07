GIW_subshed_delineation<-function(
  workspace="C:\\ScratchWorkspace\\",
  wbt_path="C:/WBT/whitebox_tools",
  dem=dem,
  depressions=depressions,
  wetland){
  
  #Set Working Directory
  setwd(paste(workspace))  
  
  #Export rater to workspace
  writeRaster(dem, paste0(workspace,"dem.tif"), overwrite=T)
  
  #Breach single-cell pits
  system(paste(paste(wbt_path), 
               "-r=BreachSingleCellPits", 
               paste0("--wd=",workspace),
               "--dem='dem.tif'", 
               "-o='dem_breachedsinglecells.tif'"))
  
  #Filter the DEM
  system(paste(paste(wbt_path), 
               "-r=EdgePreservingMeanFilter", 
               paste0("--wd=",workspace),
               "-i='dem_breachedsinglecells.tif'", 
               "-o='dem_edgepreservingfilter.tif'",
               "filter=10", 
               "threshold=100"))
  
  #Breach Analysis of the DEM
  # system(paste(paste(wbt_path), 
  #              "-r=BreachDepressions", 
  #              paste0("--wd=",workspace),
  #              "--dem=dem_edgepreservingfilter.tif", 
  #              "-o=dem_breach.tif"))
  
  #Flow Direction
  system(paste(paste(wbt_path), 
               "-r=D8Pointer", 
               paste0("--wd=",workspace),
               "--dem='dem_edgepreservingfilter.tif'", 
               "-o='fdr.tif'",
               "--out_type=sca"))
  
  #Flow Accumulation
  system(paste(paste(wbt_path), 
               "-r=DInfFlowAccumulation", 
               paste0("--wd=",workspace),
               "--dem='dem_edgepreservingfilter.tif'", 
               "-o='fac.tif'",
               "--out_type=sca"))
  
  #Read fac and fdr rasters into R environment
  fdr<-raster(paste0(workspace,"fdr.tif"))
  fac<-raster(paste0(workspace,"fac.tif"))
  
  #Extract depression of interest
  wetland_dep<-depressions
  wetland_dep[wetland_dep!=raster::extract(depressions, wetland)]<-NA
  wetland_dep<-wetland_dep*0+1
  wetland_dep[is.na(wetland_dep)]<-0
  
  #Identify pour points [all outside cells]
  p<-focal(wetland_dep, matrix(c(rep(1,4),0, rep(1,4)), nrow=3, ncol=3))
  p[p==8 | p==0]<-NA
  p<-rasterToPoints(p, spatial=T)
  p@proj4string<-dem@crs
  writeOGR(p,paste0(workspace,"."),"pp", drive="ESRI Shapefile", overwrite=T)
  
  #Define watershed pour point (max fac in depression)  
  # wetland_fac<-wetland_dep*fac
  # max_wetland_fac<-cellStats(wetland_fac, max)
  # pp<-rasterToPoints(wetland_fac, fun=function(x){x==max_wetland_fac})
  # pp<-SpatialPointsDataFrame(pp, data.frame(x=1))
  # pp@proj4string<-dem@crs
  # writeOGR(pp,paste0(workspace,"."),"pp", drive="ESRI Shapefile", overwrite=T)
  
  #Watershed Delineation
  system(paste(paste(wbt_path), 
               "-r=Watershed", 
               paste0("--wd=",workspace),
               "--d8_pntr='fdr.tif'", 
               "--pour_pts=pp.shp",
               "-o=watershed.tif"))
  
  #Import watershed shape
  watershed<-raster(paste0(workspace,"watershed.tif"))
  watershed[is.na(watershed)]<-0
  watershed<-watershed+wetland_dep
  watershed[watershed==0]<-NA
  watershed<-watershed*0+1
  
  # #Clip relevant layers to watershed
  # fac<-fac*watershed
  # fdr<-fac*watershed
  # depressions<-depressions*watershed
  # 
  # #Remove depression of interest
  # depressions[depressions==raster::extract(depressions, wetland)]<-NA
  # 
  # #Count number of depression
  # n_depression<-length(unique(depressions@data@values))-1
  # 
  # #If there are remaining depressions, remove their watershed
  # if(n_depression>0){
  #   #subset fac to reamining depressions
  #   fac<-fac*depressions
  # 
  #   #Identify boundaries
  #   r<-fac*0+1
  #   r[is.na(r)]<-0
  #   r<-focal(r, matrix(c(rep(1,4),0, rep(1,4)), nrow=3, ncol=3))
  #   r[r==8 | r==0]<-NA
  #   r<-data.frame(rasterToPoints(r))
  #   r$fac<-raster::extract(fac,SpatialPoints(r[,1:2]))
  #   r$wet<-raster::extract(depressions,SpatialPoints(r[,1:2]))
  #   
  #   #Create pour points for highest fac at each boundary
  #   r<-na.omit(r)
  #   r<- r %>% group_by(wet) %>% filter(fac==max(fac))
  #   r<- SpatialPointsDataFrame(r[,1:2], data.frame(dep=r$wet))
  #   r@proj4string<-dem@crs
  #   
  #   #Export pout points
  #   writeOGR(r,paste0(workspace,"."),"pp", drive="ESRI Shapefile", overwrite=T)
  #   
  #   #Watershed Delineation
  #   system(paste(paste(wbt_path), 
  #                "-r=Watershed", 
  #                paste0("--wd=",workspace),
  #                "--d8_pntr='fdr.tif'", 
  #                "--pour_pts=pp.shp",
  #                "-o=subshed.tif"))
  #   subshed<-raster(paste0(workspace,"subshed.tif"))
  #   
  #   #Remove subshed from watershed
  #   subshed<-subshed*0+1
  #   subshed[is.na(subshed)]<-0
  #   watershed<-watershed-subshed
  #   watershed[watershed!=1]<-NA
  # }
  
  #Export Watershed
  watershed
}
