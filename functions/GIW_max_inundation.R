GIW_max_inundation<-function(
  subshed, #wateshed raster
  dem,     #DEM for the analysis
  storage  #Storage Curve 
){
  
  #Convert to raster
  temp<-subshed*dem
  temp@crs<-dem@crs
  
  #Create Minimum Raster
  temp_min<-temp*0+minValue(temp)
  temp_min@crs<-dem@crs
  
  #Create function to return conditional raster
  Con<-function(condition, trueValue, falseValue){
    return(condition * trueValue + (!condition)*falseValue)
  }
  
  #Define relative stage for max storage capacity
  z<-min(storage$z[storage$outflow_length>0], na.omit=T)
  
  #Estimate inundated area with con function
  inundate<-Con(temp>(temp_min+z),0,1)
  inundate[inundate==0]<-NA
  
  #Export Inundate Raster
  inundate
}