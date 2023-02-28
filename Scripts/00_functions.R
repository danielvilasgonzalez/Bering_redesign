##############################
# function to plot maps R
##############################

#setwd
out_dir<-'C:/Users/Daniel.Vilas/Work/Adapting Monitoring to a Changing Seascape/slope EBS VAST/'
setwd(out_dir)

#list of sp
splist<-list.dirs('./',full.names = FALSE,recursive = FALSE)

#sp example - Gadus macrocephalus
sp<-splist[7]

#set wd
setwd(paste0('./',sp))

plot_distribution<-function(working_dir=getwd()){
  
  #required packages for plot
  library(dplyr)
  library(ggspatial)
  library(ggplot2)
  library(akgfmaps)
  library(reshape2)
  library(raster)
  library(scales)
  library(sp)
  
  # Define plot exent (through trial end error)
  # panel_extent <- data.frame(x = c(-1456559.21, -76636.05), #x = c(-1326559.21, -87636.05),
  #                            y = c(453099.5, 1894909.7))
  
  #list of models
  models<-list.dirs(getwd(),full.names = FALSE,recursive = FALSE)

  #alaska layer
  ebs_layers <- akgfmaps::get_base_layers(select.region = "ebs", set.crs = "EPSG:3338")
  ak_sppoly<-as(ebs_layers$akland, 'Spatial')

  for (m in models) {
   
    #print year to check progress
    cat(paste("\n","    ----- ", sp, " -----\n","       - ", m, " model\n"))  
    
      #load fit
      load(paste(m,'fit.RData',sep = '/'))
    
      #get dataframe
      D_gt <- drop_units(fit$Report$D_gct[,1,] )
      D_gt<-data.frame('cell'=c(1:fit$spatial_list$n_g),D_gt,check.names = FALSE)
      
      #melt df
      D_gt1<-reshape2::melt(D_gt,id=c('cell'))
      
      #make map info
      mdl <- make_map_info(Region = fit$settings$Region,
                           spatial_list = fit$spatial_list,
                           Extrapolation_List = fit$extrapolation_list)
      
      #merge map info with predictions
      D <- merge(D_gt1, mdl$PlotDF, by.x='cell', by.y='x2i')
      
      #df to spatialpoint df
      coordinates(D) <- ~ Lon + Lat
      crs(D)<-c(crs="+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
      
      #reproject coordinates for plotting purposes
      D1<-spTransform(D,'+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs')
      D2<-data.frame(D1)
      
      #plot
      p<-ggplot() +
        geom_point(data=D2, aes(Lon, Lat, color=log(as.vector(value)), group=NULL),
                   ## These settings are necessary to avoid
                   ## overlplotting which is a problem here. May need
                   ## to be tweaked further.
                   size=0.3, stroke=0,shape=16)+
        geom_polygon(data=ak_sppoly,aes(x=long,y=lat,group=group),fill = 'grey60')+
        scale_x_continuous(expand = c(0,0),
                           breaks = c(-175,-170,-165,-160))+
        scale_y_continuous(expand = c(0,0),
                           breaks = c(62,58,54))+
        #geom_polygon(data=nbs_sh,aes(x=long,y=lat,group=group),fill=NA,col='black')+
        #geom_polygon(data=ebs_sh,aes(x=long,y=lat,group=group),fill=NA,col='black')+
        coord_sf(crs = '+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs',
                 xlim = extent(D1)[1:2]-c(100000,-100000), ylim = extent(D1)[3:4]-c(100000,-100000)
                 #xlim = panel_extent$x,ylim = panel_extent$y
                 )+
        scale_color_gradientn(colours = wesanderson::wes_palette("Zissou1", 21, type = "continuous"),name=('log(kg/km²)'),
                              guide = guide_colorbar(  frame.colour = "black",ticks.colour = 'black'),na.value=rgb(1, 0, 0, 0))+
        #annotation_north_arrow(location = "tr", which_north = "true",pad_x = unit(0.01, 'points'), pad_y = unit(10, 'points'),
        #                       style = north_arrow_fancy_orienteering(line_width = 1.5, text_size =6))+
        theme(panel.border = element_rect(fill=NA),aspect.ratio = 1,panel.grid.major = element_line(color = rgb(0, 0, 0,20, maxColorValue = 285), linetype = 'dashed', size = 0.5), 
              panel.background = element_rect(fill = NA),panel.ontop = TRUE,text = element_text(size=10),
              legend.background =  element_rect(fill = "transparent", colour = "transparent"),legend.key.height= unit(12, 'points'),
              legend.key.width= unit(12,units = 'points'),axis.title = element_blank(),#legend.position = c(0.12,0.18),
              strip.background =element_rect(fill="white"),plot.title = element_text(hjust = 0.5,size=12),
              legend.title = element_text(size = 8),
              plot.margin = unit(c(-0.01,-0.01,-0.01,-0.01), "points"),legend.key = element_rect(color="black"))+
        #annotate("text", x = -256559, y = 1354909, label = "Alaska",parse=TRUE,size=7)+
        #annotate("text", x = -1176559, y = 1904909, label = "Russia",parse=TRUE,size=7)+
        labs(title=paste(sp,paste0('(',m,')') ))+
        facet_wrap(~variable)
      
      #save plot
      ggsave(paste(m,paste0(sp,"_byYear.png"),sep='/'), p, width=8, height=8)
  }
}  

#test function
plot_distribution(working_dir = getwd(),model = models)  

##############################
# plot covariate effect
##############################





########################################
# plot dev explained comparison and aic
########################################
