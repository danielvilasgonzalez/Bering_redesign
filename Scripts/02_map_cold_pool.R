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
   devtools::install_github('afsc-gap-products/coldpool')}

#load/install packages
pacman::p_load(pack_cran,character.only = TRUE)

#set working directory
mydir<-'E:/UW/Adapting Monitoring to a Changing Seascape/Resources/'
setwd(mydir)

#sbt ebs
sbt_brick<-coldpool::ebs_bottom_temperature

plot_list<-list()

#loop over years, coldpool only has data between 1985 and 2019
for (y in 1985:2019) {
  
  y<-2019
  
  #print year to check progress
  cat(paste("    ----- year", y, "-----\n"))  
  
  r<-sbt_brick[[grep(y,names(sbt_brick))]]
  
  #p<-
   gplot(r) +
    geom_tile(aes(fill=value))+
    geom_polygon(data=ak_sppoly,aes(x=long,y=lat,group=group),fill = 'grey60')+
    #geom_polygon(data=eez_sh22,aes(x=long,y=lat,group=group),fill=NA,color='purple')+
    #geom_polygon(data=eez_sh11,aes(x=long,y=lat,group=group),fill=NA,color='purple')+
    #geom_polygon(data=eez_sh33,aes(x=long,y=lat,group=group),fill=NA,color='grey40')+
    scale_x_continuous(expand = c(0,0),position = 'bottom',
                       breaks = c(-175,-170,-165,-160,-155),sec.axis = dup_axis())+
    #geom_polygon(data=nbs_sh,aes(x=long,y=lat,group=group),fill=NA,col='black')+
    geom_polygon(data=ebs_sh,aes(x=long,y=lat,group=group),fill=NA,col='black')+
    coord_sf(crs = '+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs',
             xlim = panel_extent$x,
             ylim = c(403099.5, 1754909.7),
             label_axes = "-NE-")+
    scale_fill_binned(limits=c(-2,2),na.value = 'white',low = "#007c9b", high = "#c1f2fe",name='SBT',guide = guide_colorbar(frame.colour = 'black',ticks.colour = NA))+
    # scale_fill_gradient2(limits=c(-2,2),na.value = 'white',low = '#03254c',
    #                                           mid='#187bcd',
    #                                            high = '#d0efff',guide = guide_colorbar(frame.colour = 'black',ticks.colour = 'black'))+
    # guides(fill = guide_legend(direction = "horizontal",
    #                            title.position = "top",
    #                            label.position = "bottom",
    #                            label.hjust = 0.5,
    #                            label.vjust = 1),title='SBT')+
    #scale_shape_manual(values=c('98'=1,'143'=2))+
    # scale_fill_gradient2(low = '#c1f2fe',
    #                      #mid='#5fdfff',
  #                      high = '#007c9b',
  #                      limits=c(0,200),
  #                      na.value = 'white',
  #                      name='depth (m)',
  #                      guide = guide_colorbar(frame.colour = 'black',ticks.colour = 'black'))+
  # #annotation_north_arrow(location = "tr", which_north = "true",pad_x = unit(0.01, 'points'), pad_y = unit(5, 'points'),
  #                       style = north_arrow_fancy_orienteering(line_width = 1.5, text_size =6))+
  theme(aspect.ratio = 1,panel.grid.major = element_line(color = rgb(0, 0, 0,20, maxColorValue = 285), linetype = 'dashed', size = 0.5),
        panel.background = element_rect(fill = NA),panel.ontop = TRUE,text = element_text(size=10),
        legend.background =  element_rect(fill = "transparent", colour = "transparent"),legend.key.height= unit(20, 'points'),
        legend.key.width= unit(20, 'points'),axis.title = element_blank(),legend.position = c(0.12,0.50),
        panel.border = element_rect(fill = NA, colour = 'black'),legend.key = element_rect(color="black"),
        axis.text = element_text(color='black'),legend.spacing.y = unit(10, 'points'))+
    theme(    axis.text.y.right = element_text(hjust= 0.1 ,margin = margin(0,0,0,-25, unit = 'points'),color='black'),
              axis.text.x = element_text(vjust = 6, margin = margin(-7,0,0,0, unit = 'points'),color='black'),
              axis.ticks.length = unit(-5,
                                       "points"))+
    #annotate("text", x = -256559, y = 1354909, label = "Alaska",parse=TRUE,size=7)+
    scale_y_continuous(expand = c(0,0),position = 'right',sec.axis = dup_axis())+
    #annotate("text", x = -1376559, y = 814900, label = "italic('Eastern')",parse=TRUE,size=9)+
    #annotate("text", x = -1376559, y = 694900, label = "italic('Bering Sea')",parse=TRUE,size=9)
    annotate("text", x = -556559, y = 1624909, label = y,size=9)
  
  plot_list[[paste0(y)]]<-p
  
  
  
}


tiff('./Daniel Bering Sea/coldpool_1.tiff',height=2500,width = 2500,res = 250)
cowplot::plot_grid(plotlist = plot_list[1:9],nrow = 3,ncol = 3)
dev.off()

tiff('./Daniel Bering Sea/coldpool_2.tiff',height=2500,width = 2500,res = 250)
cowplot::plot_grid(plotlist = plot_list[10:18],nrow = 3,ncol = 3)
dev.off()

tiff('./Daniel Bering Sea/coldpool_3.tiff',height=2500,width = 2500,res = 250)
cowplot::plot_grid(plotlist = plot_list[19:27],nrow = 3,ncol = 3)
dev.off()

tiff('./Daniel Bering Sea/coldpool_4.tiff',height=2500,width = 2500,res = 250)
cowplot::plot_grid(plotlist = plot_list[28:35],nrow = 3,ncol = 3)
dev.off()
