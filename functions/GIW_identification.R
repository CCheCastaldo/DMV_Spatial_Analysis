GIW_identification<-function(
  dem=dem, 
  min_size=100, #number of cells (for now)
  iterations=100, #number of Monte Carlo iterations
  dem_rmse=0.0607, #RMSE of 18.5 cm
  workspace="C:\\ScratchWorkspace\\", 
  wbt_path="C:/WBT/whitebox_tools"){
  
  #Set working directory
  setwd(paste(workspace))
  
  #Write raster to ScratchWorkspace
  dem<-na.omit(dem)
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
  
  #Identify depressions using WhiteBox GAT Monte Carlo Approach
  system(paste(paste(wbt_path),
               "-r=StochasticDepressionAnalysis", 
               paste0("--wd=",workspace),
               "--dem='dem_edgepreservingfilter.tif'", 
               "-o='depression.tif'",
               paste0("--rmse=",dem_rmse),
               paste0("--iterations=",iterations)))
  
  #Reclass raster
  system(paste(paste(wbt_path), 
               "-r=Reclass", 
               paste0("--wd=",workspace),
               "-i='depression.tif'", 
               "-o='reclass.tif'",
               "--reclass_vals='0;0;0.80';1;0.80;1"))
  
  #Identify clusters of inundated cells
  system(paste(paste(wbt_path), 
               "-r=Clump",
               paste0("--wd=",workspace),
               "-i='reclass.tif'", 
               "-o='group.tif'",
               "--diag", 
               "--zero_back"))
  
  #Identify Clusters greater than 100m2
  r<-raster(paste0(workspace,"group.tif")) #Read raster
  r_pnts<-data.frame(rasterToPoints(r))        #Convert to XYZ pnts
  r_remove<-r_pnts %>% group_by(group) %>% tally() #sum # of raster cells for each group
  r_remove<-r_remove$group[r_remove$n>min_size]
  r_pnts$group[r_pnts$group %in% r_remove]<-0
  r_pnts<-r_pnts[r_pnts$group!=0,]
  r_pnts<-SpatialPoints(r_pnts[,1:2])
  r_remove<-rasterize(r_pnts, r, 0)
  r_remove[is.na(r_remove)]<-1
  r<-r*r_remove
  r[r==0]<-NA
  
  #Export depression raster
  r
}
