####################################################################
####################################################################
##
##    Create a cold pool extension map
##    from a) coldpool package and b) ROMS Bering Sea
##    Daniel Vilas (danielvilasgonzalez@gmail.com/dvilasg@uw.edu)
##
####################################################################
####################################################################

#clear all objects
rm(list = ls(all.names = TRUE)) 
#free up memrory and report the memory usage
gc() 

#libraries from cran to call or install/load
pack_cran<-c('rasterVis','scales','rnaturalearth','cowplot','ggplot2')

#install pacman to use p_load function - call library and if not installed, then install
if (!('pacman' %in% installed.packages())) {
  install.packages("pacman")}

#install coldpool to extract SBT for the EBS
if (!('coldpool' %in% installed.packages())) {
   devtools::install_github('afsc-gap-products/coldpool')};library(coldpool)

#install akgfmaps to extract shapefile of Alaska
if (!('akgfmaps' %in% installed.packages())) {
  devtools::install_github('afsc-gap-products/akgfmaps')};library(akgfmaps)

#load/install packages
pacman::p_load(pack_cran,character.only = TRUE)

#set working directory
mydir<-'E:/UW/Adapting Monitoring to a Changing Seascape/Resources/'
setwd(mydir)

# Define plot extent (through trial end error) units km
panel_extent <- data.frame(x = c(-1716559.21, -77636.05), #x = c(-1326559.21, -87636.05),
                           y = c(483099.5, 2194909.7)) #y = c(533099.5, 1894909.7))

#Alaska layer
ebs_layers <- akgfmaps::get_base_layers(select.region = "ebs", set.crs = "EPSG:3338")
ak_sppoly<-as(ebs_layers$akland, 'Spatial')

#shapefile EBS
ebs_sh<-rgdal::readOGR(dsn='./Shapefiles/',layer = 'EBSshelfThorson')

#sbt ebs
sbt_brick<-coldpool::ebs_bottom_temperature

#names
years<-sort(as.numeric(gsub(".*?([0-9]+).*", "\\1", names(sbt_brick))))

#list to store values
plot_list<-list()

#loop over years, coldpool only has data between 1982 and 2022
for (y in min(years):max(years)) {

  #y<-1982
  #print year to check progress
  cat(paste("    ----- year", y, "-----\n"))  
  
  #2020 empty, so all NAs
  if (y == 2020) {
    r<-sbt_brick[[grep(y-1,names(sbt_brick))]]
    r[r]<-NA
  } else {
    r<-sbt_brick[[grep(y,names(sbt_brick))]]
  }
  
  #plot
  p<-
   gplot(r) +
    geom_tile(aes(fill=value))+
    geom_polygon(data=ak_sppoly,aes(x=long,y=lat,group=group),fill = 'grey60')+
    scale_x_continuous(expand = c(0,0),position = 'bottom',
                       breaks = c(-175,-170,-165,-160,-155),sec.axis = dup_axis())+
    geom_polygon(data=ebs_sh,aes(x=long,y=lat,group=group),fill=NA,col='black')+
    coord_sf(crs = '+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs',
             xlim =c(-1716559.21,57636.05),
             ylim = c(403099.5, 1754909.7),
             label_axes = "-NE-")+
    scale_fill_binned(limits=c(-2,2),na.value = 'white',low = "#007c9b", high = "#c1f2fe",name='SBT (°C)',guide = guide_colorbar(frame.colour = 'black',ticks.colour = NA))+
    theme(aspect.ratio = 1,panel.grid.major = element_line(color = rgb(0, 0, 0,20, maxColorValue = 285), linetype = 'dashed', linewidth = 0.5),
        panel.background = element_rect(fill = NA),panel.ontop = TRUE,text = element_text(size=10),
        legend.background =  element_rect(fill = "transparent", colour = "transparent"),legend.key.height= unit(10, 'points'),
        legend.key.width= unit(10, 'points'),axis.title = element_blank(),legend.position = c(0.15,0.35),
        panel.border = element_rect(fill = NA, colour = 'black'),legend.key = element_rect(color="black"),
        axis.text = element_text(color='black'),legend.spacing.y = unit(10, 'points'),
        axis.text.y.right = element_text(hjust= 0.1 ,margin = margin(0,0,0,-25, unit = 'points'),color='black'),
        axis.text.x = element_text(vjust = 6, margin = margin(-7,0,0,0, unit = 'points'),color='black'),
        axis.ticks.length = unit(-5,"points"))+
    scale_y_continuous(expand = c(0,0),position = 'right',sec.axis = dup_axis())+
    annotate("text", x = -1406559, y = 1624909, label = y,size=8)

  #remove legend from plots except first (topleft) in the panel
  if (!(y %in% years[c(1,17,33)])) {
    p<-p+theme(legend.position = 'none')}
  
  #store plot on list
  plot_list[[paste0(y)]]<-p
  
  
  
}

#save plot 16 first years (4x4)
tiff(paste0('E:/UW/Adapting Monitoring to a Changing Seascape/Figures/coldpool_',years[1],'_',years[16],'.tiff'),height=2500,width = 2500,res = 250)
cowplot::plot_grid(plotlist = plot_list[1:16],nrow = 4,ncol = 4)
dev.off()

#save plot 16 next years (4x4)
tiff(paste0('E:/UW/Adapting Monitoring to a Changing Seascape/Figures/coldpool_',years[17],'_',years[33],'.tiff'),height=2500,width = 2500,res = 250)
cowplot::plot_grid(plotlist = plot_list[17:33],nrow = 4,ncol = 4)
dev.off()

#save plot last years (4x4)
tiff(paste0('E:/UW/Adapting Monitoring to a Changing Seascape/Figures/coldpool_',years[34],'_',years[length(years)],'.tiff'),height=2500,width = 2500,res = 250)
cowplot::plot_grid(plotlist = plot_list[34:length(years)],nrow = 4,ncol = 4)
dev.off()



