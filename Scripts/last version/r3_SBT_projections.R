####################################################################
####################################################################
##
##    check temperature in EBS and NBS
##    create raster of SBT projections
##    Daniel Vilas (danielvilasgonzalez@gmail.com/dvilasg@uw.edu)
##
####################################################################
####################################################################

#clear all objects
rm(list = ls(all.names = TRUE)) 
#free up memrory and report the memory usage
gc() 

#libraries from cran to call or install/load
pack_cran<-c('ncdf4','raster','FNN','lubridate','rgeos','scales','rnaturalearth','grid','ggplot2','rasterVis','ggthemes','reshape2','ragg',"wesanderson")

#install pacman to use p_load function - call library and if not installed, then install
if (!('pacman' %in% installed.packages())) {
  install.packages("pacman")}

#load/install packages
pacman::p_load(pack_cran,character.only = TRUE)

#setwd - depends on computer using
#out_dir<-'C:/Users/Daniel.Vilas/Work/Adapting Monitoring to a Changing Seascape/' #NOAA laptop  
out_dir<-'/Users/daniel/Work/Adapting Monitoring to a Changing Seascape/' #mac
#out_dir<-'/Users/daniel/Work/VM' #VM
setwd(out_dir)

#range years of data
sta_y<-1982
end_y<-2022

# #get files from google drive and set up
# files<-googledrive::drive_find()
# 32 #for dvilasg@uw.edu
# 
# #get id shared folder from google drive
# id.bering.folder<-files[which(files$name=='Bering redesign RWP project'),'id']

#selected species
spp<-c('Limanda aspera',
       'Gadus chalcogrammus',
       'Gadus macrocephalus',
       'Atheresthes stomias',
       'Reinhardtius hippoglossoides',
       'Lepidopsetta polyxystra',
       'Hippoglossoides elassodon',
       'Pleuronectes quadrituberculatus',
       'Hippoglossoides robustus',
       'Boreogadus saida',
       'Eleginus gracilis',
       'Anoplopoma fimbria',
       'Chionoecetes opilio',
       'Paralithodes platypus',
       'Paralithodes camtschaticus')

# #list of files and folder
# files.1<-googledrive::drive_ls(id.bering.folder$id)
# id.data<-files.1[which(files.1$name=='data raw'),'id']
# files.2<-googledrive::drive_ls(id.data$id)

#check temperature example for Pcod
#for (sp in spp) {
  
  #sp
  sp<-'Gadus macrocephalus'
  
  #open data
  df1<-readRDS(paste0('./data processed/species/',sp,'/data_geostat_temp.rds'))
  df1<-subset(df1,year %in% 1982:2022)
  
  #select rows and rename
  df2<-df1[,c("lat_start","lon_start","year",'scientific_name','weight_kg','effort','depth_m','LogDepth',"ScaleLogDepth",'Scalebottom_temp_c','bottom_temp_c','survey_name')]
  colnames(df2)<-c('Lat','Lon','Year','Species','CPUE_kg','Effort','Depth','LogDepth','ScaleLogDepth','ScaleBotTemp','BotTemp','Region')
  
  #data geostat for EBS and NBS
  df3<-subset(df2,Region %in% c("Eastern Bering Sea Crab/Groundfish Bottom Trawl Survey",
                                "Northern Bering Sea Crab/Groundfish Survey - Eastern Bering Sea Shelf Survey Extension"))
  #range og years
  yrs_region<-range(df3$Year)
  
  ######################################
  # SBT ANALYSIS
  ######################################
  
  #mean SBT by year
  yearagg.df <- aggregate(data = df3, BotTemp ~ Year, mean)
  yearagg.df$TempAnomaly<-NA
  
  #calculate SBT anomaly
  for (i in 2:nrow(yearagg.df)) {
    yearagg.df[i,'TempAnomaly']<-yearagg.df[i,'BotTemp']-yearagg.df[i-1,'BotTemp']
  }
  
  #plot 1
  print(
    ggplot() +
      geom_line(data=yearagg.df, aes(x=Year, y=BotTemp),linetype='dashed')+
      geom_point(data=yearagg.df, aes(x=Year, y=BotTemp,color=BotTemp),size=2)+
      geom_bar(data=yearagg.df, aes(x=Year, y=TempAnomaly,fill=TempAnomaly),color='black',stat="identity",position = position_dodge(0.9))+
      scale_colour_gradient2(low = 'darkblue',high='darkred',midpoint = 2.5)+
      scale_fill_gradient2(low = 'darkblue',high='darkred',midpoint = 0)+
      #xlab(label = 1982:2022)+
      scale_x_continuous(breaks=c(1982:2022))+
      theme_bw()+
      labs(title=paste0(sp))+
      theme(panel.grid.minor.y = element_blank())
  )
  
  #plot 2
  print(
    ggplot() +
      geom_density(data=df3, aes(x=BotTemp, group=Year))+
      geom_vline(data=yearagg.df,aes(xintercept=BotTemp,group=Year, color=BotTemp),
                   linetype="dashed", size=1)+
      scale_colour_gradient2(low = 'darkblue',high='darkred',midpoint = 2.5)+
      geom_text(data=yearagg.df,aes(label=paste0('mean BotTemp = ',round(BotTemp,digits = 2)),x = 10,y=0.43))+
      facet_wrap(~Year)+
      labs(title=paste0(sp))+
      theme_bw())
#}

######################################
# SBT PROJECTIONS
######################################

#create directory
dir.create('./tables/',showWarnings = FALSE)
  
#create SBT projections
df_sbt<-data.frame(y2022=c(0,0,0,0,0,0,0,0,0,0,0,0),
                   y2023=c(0,-0.5,0.5,1,1,-1,2,-2,3,-3,4,5),
                   y2024=c(0,-1,1,2,0,0,0,0,0,0,0,0),
                   y2025=c(0,-1.5,1.5,3,1,-1,2,-2,3,-3,4,5),
                   y2026=c(0,-2,2,4,0,0,0,0,0,0,0,0),
                   y2027=c(0,-2.5,2.5,5,1,-1,2,-2,3,-3,4,5),
                   Scenario=c('status quo','gradually cold','gradually warm','severe gradually warm',
                              'warm low variation','cold low variation','warm moderate variation',
                              'cold moderate variation','warm high variation','cold high variation',
                              'warm very high variation','warm extreme variation'))
  
#remove cold SBT scenarios
df_sbt<-df_sbt[-grep('cold',df_sbt$Scenario),]    
df_sbt$sbt_n<-1:nrow(df_sbt)

#save SBT table
save(df_sbt,file = './tables/SBT_projection.RData')

#load grid
load('./data processed/grid_EBS_NBS.RData')

#subset to year 2022
df2022<-subset(grid.ebs_year,Year==2022)

#create spatial object from df
coordinates(df2022) <- ~ Lon + Lat
proj4string(df2022) <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0") 
df2022<-spTransform(df2022,CRSobj = CRS('+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs'))

#n cells
xycells<-as.integer(sqrt(dim(df2022)[1]))

# create a template raster
r1 <- raster(ext=extent(df2022),nrow=xycells,ncol=xycells)

#create raster
r2<-rasterize(df2022, r1 ,field='Temp')
crs(r2) <- CRS('+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs')
#plot(r1)

#color palette
pal<-wesanderson::wes_palette('Zissou1',21,type='continuous')

#year
yrs<-2023:2027

#dir.create
dir.create('./data processed/SBT projections/')

#create folder
dir.create('./figures/SBT/')

#Alaska layer
ebs_layers <- akgfmaps::get_base_layers(select.region = "ebs", set.crs = "EPSG:3338")
ak_sppoly<-as(ebs_layers$akland, 'Spatial')

# Define plot extent (through trial end error) units km
panel_extent <- data.frame(x = c(-1516559.21, -77636.05), #x = c(-1326559.21, -87636.05),
                           y = c(483099.5, 1894909.7)) #y = c(533099.5, 1894909.7))

#loop over scenarios and projected years
for (sbt in unique(df_sbt$sbt_n)) {

  #sbt<-8
  
  plot_list<-list()
  
  colour_breaks <- c(1,1.5, 2,2.5, 9)
  #colours <- c('darkblue',"lightblue", "grey90", "red",'darkred')
  
  #p<-
    gplot(r2)+
    geom_tile(aes(fill=value))+
    geom_polygon(data=ak_sppoly,aes(x=long,y=lat,group=group),color='black',linewidth=0.2,fill = 'grey80')+
    #geom_polygon(data=eez_sh33,aes(x=long,y=lat,group=group),fill=NA,color='grey40')+
    scale_x_continuous(expand = c(0,0),position = 'bottom',
                       breaks = c(-175,-170,-165,-160,-155),sec.axis = dup_axis())+
    coord_sf(crs = '+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs',
             xlim = panel_extent$x,
             ylim = panel_extent$y,
             #lim = c(panel_extent$x[1]+200000,panel_extent$x[2]),
             #ylim = c(panel_extent$y[1],panel_extent$y[2]-200000),
             label_axes = "-NE-")+
      scale_fill_gradientn(
        limits  = c(-1,15),
        colours = colours[c(1, seq_along(colours), length(colours))],
        values  = c(0, scales::rescale(colour_breaks, from = c(-1,15)), 1))+
    # scale_fill_gradient2(low = 'blue',midpoint = 4,#c1f2fe',
    #                      high = 'red',#'#007c9b',
    #                      limits=c(-1,15),#oob = scales::squish,breaks=c(0,50,100,200),
    #                      #labels=c('0','50','100',paste0('200 - ',round(maxValue(ak_bathy_4)))),
    #                      na.value = 'white',
    #                      name='SBT',
    #                      guide = guide_colorbar(frame.colour = 'black',ticks.colour = 'black'))+
    theme(aspect.ratio = 1,panel.grid.major = element_line(color = rgb(0, 0, 0,20, maxColorValue = 285), linetype = 'dashed', linewidth =  0.5),
          panel.background = element_rect(fill = NA),panel.ontop = TRUE,text = element_text(size=10),
          legend.background =  element_rect(fill = "transparent", colour = "transparent"),legend.key.height= unit(25, 'points'),
          legend.key.width= unit(25, 'points'),axis.title = element_blank(),legend.position = 'none', #c(0.12,0.47)
          panel.border = element_rect(fill = NA, colour = 'black'),legend.key = element_rect(color="black"),
          axis.text = element_text(color='black'),legend.spacing.y = unit(10, 'points'),
          axis.text.y.right = element_text(hjust= 0.1 ,margin = margin(0,7,0,-25, unit = 'points'),color='black'),
          axis.text.x = element_text(vjust = 6, margin = margin(-7,0,7,0, unit = 'points'),color='black'),
          axis.ticks.length = unit(-5,"points"),plot.title = element_text(hjust = 0.01, vjust = -7, size = 14),
          plot.margin = margin(0.01,0.01,0.01,0.01))+
    
    labs(title = '2022')+
    scale_y_continuous(expand = c(0,0),position = 'right',sec.axis = dup_axis())
  
  plot_list[[1]]<-p
    
  legend1 <- cowplot::get_legend( 
    p + 
      theme(legend.position = "bottom")+ guides(fill=guide_colorbar(title.position='top', title.hjust = 0.5,ticks.colour = 'black',frame.colour = 'black'))
  ) 
  
  #plot(legend1)
  
  #create stack
  sbt_stack<-stack()
  
  #add layer to stack
  sbt_stack<-addLayer(sbt_stack,r2)
  names(sbt_stack)<-'y2022'
  
  #loop over projected years
  for (y in 1:5) {
    #y<-1
    
    #print species to check progress
    cat(paste(" #############  Scenario", sbt, "year",y, " #############\n"))
    
    #modify raster y2022
    r3<-r2
    names(r3)<-paste0('y',yrs[y])
    values(r3)<-values(r3)+df_sbt[which(df_sbt$sbt_n==sbt),paste0('y',yrs[y])]
    
    #add layer to stack
    sbt_stack<-addLayer(sbt_stack,r3)
    #r4<-as.data.frame(r3,xy=TRUE)
    #colnames(r4)[3]<-'value'
    #r4<-r4[complete.cases(r4$value),] 
    
    #ggplot()+
      #geom_raster(data=r4,aes(x=x,y=y,fill=value))+
    p<-
    gplot(r3)+
      geom_tile(aes(fill=value))+
      geom_polygon(data=ak_sppoly,aes(x=long,y=lat,group=group),color='black',linewidth=0.2,fill = 'grey80')+
      #geom_polygon(data=eez_sh33,aes(x=long,y=lat,group=group),fill=NA,color='grey40')+
      scale_x_continuous(expand = c(0,0),position = 'bottom',
                         breaks = c(-175,-170,-165,-160,-155),sec.axis = dup_axis())+
      coord_sf(crs = '+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs',
               xlim = panel_extent$x,
               ylim = panel_extent$y,
               #lim = c(panel_extent$x[1]+200000,panel_extent$x[2]),
               #ylim = c(panel_extent$y[1],panel_extent$y[2]-200000),
               label_axes = "-NE-")+
      scale_fill_gradient2(low = 'blue',midpoint = 2,#c1f2fe',
                           high = 'red',#'#007c9b',
                           limits=c(-1,15),#oob = scales::squish,breaks=c(0,50,100,200),
                           #labels=c('0','50','100',paste0('200 - ',round(maxValue(ak_bathy_4)))),
                           na.value = 'white',
                           name='SBT',
                           guide = guide_colorbar(frame.colour = 'black',ticks.colour = 'black'))+
      theme(aspect.ratio = 1,panel.grid.major = element_line(color = rgb(0, 0, 0,20, maxColorValue = 285), linetype = 'dashed', linewidth =  0.5),
            panel.background = element_rect(fill = NA),panel.ontop = TRUE,text = element_text(size=10),
            legend.background =  element_rect(fill = "transparent", colour = "transparent"),legend.key.height= unit(25, 'points'),
            legend.key.width= unit(25, 'points'),axis.title = element_blank(),legend.position = 'none', #c(0.12,0.47)
            panel.border = element_rect(fill = NA, colour = 'black'),legend.key = element_rect(color="black"),
            axis.text = element_text(color='black'),legend.spacing.y = unit(10, 'points'),
            plot.margin = margin(0.01,0.01,0.01,0.01), 
            axis.text.y.right = element_text(hjust= 0.1 ,margin = margin(0,7,0,-25, unit = 'points'),color='black'),
            axis.text.x = element_text(vjust = 6, margin = margin(-7,0,7,0, unit = 'points'),color='black'),
            axis.ticks.length = unit(-5,"points"),plot.title = element_text(hjust = 0.01, vjust = -7, size = 14))+
      labs(title = yrs[y])+
      scale_y_continuous(expand = c(0,0),position = 'right',sec.axis = dup_axis())
      
    plot_list[[y+1]]<-p
    
    
  }

  
  pgrid1<-cowplot::plot_grid(plotlist = plot_list, nrow = 1)
  #pgrid2<-cowplot::plot_grid(plotlist = plot_list_n, nrow = 2)
  title <- cowplot::ggdraw() + cowplot::draw_label(paste0(df_sbt[sbt,'Scenario'],' (SBT',sbt,')'), fontface='bold')
  
  #save plots
  ragg::agg_png(paste0('./figures/sbt',sbt,'.png'), width = 20, height = 5, units = "in", res = 300)
  print(cowplot::plot_grid(title,pgrid1, legend1, nrow = 3, rel_heights = c(.05,1, .3)))
  dev.off()
  
  
  
  
  #sbt remove 2022
  sbt_stack<-sbt_stack[[paste0('y',2023:2027)]]
  
  
  
  
  # #save raster
  # writeRaster(sbt_stack, paste0('data processed/SBT projections/SBT_',sbt,'_',paste0(range(yrs)[1],'-',range(yrs)[2]),'.grd'),overwrite=TRUE)
  # 
  # #save plot
  # ragg::agg_png(paste0('./figures/SBT/SBT_',sbt,'.png'), width = 20, height = 4.6, units = "in", res = 300)
  # print(
  #   levelplot(sbt_stack,
  #             layout=c(6, 1),
  #             col.regions=pal,
  #             main=paste('Scenario',
  #                        df_sbt[which(df_sbt$sbt_n==sbt),"Scenario"]),
  #             contour=TRUE))
  # dev.off()
}

#projection name
df_sbt$Scenario<-paste0('SBT',1:nrow(df_sbt),'_',df_sbt$Scenario)

#reshape df
df_sbt1<-reshape2::melt(df_sbt[,1:7],id='Scenario')
df_sbt1$Scenario <- factor(df_sbt1$Scenario, levels = df_sbt$Scenario)

#plot
p<-
  ggplot()+
  geom_line(data=df_sbt1,aes(x=variable,y=value,group=Scenario,color=Scenario))+
  scale_color_tableau(palette = 'Tableau 10')+
  #scale_color_manual(values=c("#ea5545", "#f46a9b", "#ef9b20", "#edbf33", "#ede15b", "#bdcf32", "#87bc45", "#27aeef", "#b33dc6"))+
  theme_bw()+
  scale_x_discrete(labels=c(2022:2027))+
  labs(color='SBT scenarios',x='',y='SBT anomaly (°C)')

ragg::agg_png(paste0('./figures/SBT_projections.png'), width = 7, height = 4, units = "in", res = 300)
print(p)
dev.off()


# #Table  
# colnames(df_sbt)<-c(2022:2027,'Scenario')
# df_sbt<-df_sbt[order(df_sbt$Scenario),]
# table<-ggpubr::ggtexttable(df_sbt, rows = NULL)

