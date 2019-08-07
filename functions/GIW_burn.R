#Define variables (for eventual function)
GIW_burn<-function(
  dem,   #DEM in question 
  burn){ #Wetland Polygon
  
  #Create DEM mask
  mask<-rasterize(burn, dem, 1)
  dem_mask<-mask*dem
  
  #Create minimum raster
  dem_min<-cellStats(dem_mask, min, na.rm=T)
  dem_min<-dem_mask*0+dem_min
  dem_min[is.na(dem_min)]<-0
  
  #Replace masked location with min raster value
  dem_mask<-dem_mask*0
  dem_mask[is.na(dem_mask)]<-1
  dem<-dem*dem_mask+dem_min
  
  #Export DEM
  dem
}