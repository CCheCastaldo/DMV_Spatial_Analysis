GIW_stage_storage<-function(
  subshed, #wateshed raster
  dem,     #DEM for the analysis
  z_max,   #Max inundation depth in map units
  dz       #Inundation increments in map units
){
  
  #Convert to raster
  temp<-subshed*dem
  temp@crs<-dem@crs
  
  #Create Minimum Raster
  temp_min<-temp*0+minValue(temp)
  temp_min@crs<-dem@crs
  
  #Create minimum point
  pnt_min<-cellStats(temp, min)
  pnt_min<-rasterToPoints(temp,spatial=T, fun=function(x) x ==pnt_min)

  #Create function to return conditional raster
  Con<-function(condition, trueValue, falseValue){
    return(condition * trueValue + (!condition)*falseValue)
  }
  
  #Create function to calcluate inundation area, volume, and spill boundary length
  inundate<-function(z){
    
    #define inundated area [area connected to pnt_min]
    area<-Con(temp>(temp_min+z),0,1)
    area<-clump(area)
    group<-raster::extract(x=area, y=pnt_min)
    area[area!=group]<-0
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
  
  #print dataframe
  df
}