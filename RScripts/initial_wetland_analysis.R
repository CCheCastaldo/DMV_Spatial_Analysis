#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Title: Initial Wetand Analysis 
#Coder: C. Nathan Jones (cnjones7@ua.edu)
#Date: 10/2/2019
#Purpose: Examine hydrogeomorphic features accross Palmer Lab wetland sites
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#1.0 Setup workspace============================================================
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Clear memory
rm(list=ls(all=TRUE))

#Defin relevant working directories
data_dir<-"//nfs/palmer-group-data/choptank/Nate/Spatial_Analysis/"

#Download packages 
library(leaflet)
library(sf)
library(raster)
library(corrplot)
library(tidyverse)

#read data from spatial analysis
giws<-st_read(paste0(data_dir, "III_Products/giws.shp"))
colnames(giws)<-c("WetID" ,"z","merge_to","spawned","root_giw","flow_to",
                  "subshed_area_m2","watershed_area_m2","volume_m3",
                  "wetland_hsc_cm","merged_hsc_cm","watershed_hsc_cm","a_axis_length_m",
                  "perimeter_m","area_m2","p_a_ratio","hand_m",
                  "mean_elevation_m","hans_m","wet_order", "geometry") 
p<-st_crs(giws)

#define wetland sites
sites<-read_csv(paste0(data_dir, "I_Data/wetland_locations.csv"))
sites<-st_as_sf(sites, coords = c("POINT_X", "POINT_Y"), crs = 2248)
sites<-st_transform(sites, crs = p)

#define intersections
df<- st_join(sites, giws) %>% filter(spawned==0) %>% filter(wet_order<10) %>% filter(Wetland!="DF") %>% st_drop_geometry()
  
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#2.0 Initial Analysis ==========================================================
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Plot all the things------------------------------------------------------------
df<-df %>% select(c(
                "subshed_area_m2","watershed_area_m2","volume_m3",
                "wetland_hsc_cm","watershed_hsc_cm","a_axis_length_m",
                "perimeter_m","area_m2","p_a_ratio","hand_m",
                "mean_elevation_m","hans_m","wet_order")) 
m<-cor(df)
corrplot(m, method="circle", type="lower")


cor.mtest <- function(mat, ...) {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat<- matrix(NA, n, n)
  diag(p.mat) <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(mat[, i], mat[, j], ...)
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
    }
  }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  p.mat
}

p.mat <- cor.mtest(df)
corrplot(m, type="upper", order="hclust", 
         p.mat = p.mat, sig.level = 0.01)

#Plot---------------------------------------------------------------------------
plot(df$wet_order, df$watershed_area_m2, log="y", xlab="Wetland Order", ylab='Watershed Area [m2]', pch=19, col="dark blue", cex=2, cex.lab=14/12, cex.axis=10/12)
plot(df$wet_order, df$volume_m3, log="y", xlab="Wetland Order", ylab='Wetland Storage [m3]', pch=19, cex=2, col="dark blue", cex.lab=14/12, cex.axis=10/12)
plot(df$wet_order, df$wetland_hsc_cm, log="y", xlab="Wetland Order", ylab='Wetland Storage [cm]', pch=19, cex=2,col="dark blue", cex.lab=14/12, cex.axis=10/12)
plot(df$wet_order, df$watershed_hsc_cm, log="y", xlab="Wetland Order", ylab='Watershed Storage [cm]', pch=19, cex=2,col="dark blue", cex.lab=14/12, cex.axis=10/12)
plot(df$wet_order, df$hand_m, xlab="Wetland Order", ylab='HAND [m]', pch=19, col="dark blue",  cex=2,cex.lab=14/12, cex.axis=10/12)
plot(df$wet_order, df$hans_m, xlab="Wetland Order", ylab='HANS [m]', pch=19, col="dark blue",  cex=2,cex.lab=14/12, cex.axis=10/12)

