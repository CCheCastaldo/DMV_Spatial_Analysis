#Quick elevation profile for agu!

#Clear memory
rm(list=ls(all=TRUE))

#Add correct libraries
library('tidyverse')
library('raster')
library('sf')

#Defin relevant working directories
data_dir<-"//nfs/palmer-group-data/choptank/Nate/Group_Maps/201805_SWS_MP/SpatialData/"
list.files(data_dir)

#define dem
dem<-raster(paste0(data_dir, 'dem_clip'))
dem<-rasterToPoints(dem)

#Grab dem data
df<-dem %>% as_tibble(.) %>% sample_n(10000)

#plot in ggplot
df %>% 
  ggplot(aes(x=dem_clip)) + 
    geom_histogram(bins=100)+
    coord_flip()+
    theme_bw()+
      labs(x = "Elevation [m]", 
           y = NULL) +
      theme(axis.text=element_text(size=14), 
            axis.title=element_text(size=18), 
            axis.title.x=element_blank(),
            axis.text.x=element_blank(),
            axis.ticks.x=element_blank())

ggsave(paste0("//nfs/njones-data/agu.png"), width=5, height=4, units="in")

