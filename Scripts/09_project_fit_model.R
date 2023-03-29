####################################################################
####################################################################
##
##    Project model example
##    Daniel Vilas (danielvilasgonzalez@gmail.com/dvilasg@uw.edu)
##
####################################################################
####################################################################

#clear all objects
rm(list = ls(all.names = TRUE)) 
#free up memrory and report the memory usage
gc() 

#libraries from cran to call or install/load
pack_cran<-c("splines")

#install pacman to use p_load function - call library and if not installed, then install
if (!('pacman' %in% installed.packages())) {
  install.packages("pacman")}

#install VAST if necessary
if (!('VAST' %in% installed.packages())) {
  devtools::install_github("james-thorson/VAST@main", INSTALL_opts="--no-staged-install")};library(VAST)

#install akgfmaps if necessary
if (!('akgfmaps' %in% installed.packages())) {
  devtools::install_github("afsc-gap-products/akgfmaps", build_vignettes = TRUE)};library(akgfmaps)

#load/install packages
pacman::p_load(pack_cran,character.only = TRUE)

#setwd
out_dir<-'C:/Users/Daniel.Vilas/Work/Adapting Monitoring to a Changing Seascape/'
setwd(out_dir)

#version VAST (cpp)
version<-'VAST_v13_1_0'

##############################
#OBJECTS TO PLOT
#############################
# Define plot extent (through trial end error) units km
panel_extent <- data.frame(x = c(-1716559.21, -77636.05), #x = c(-1326559.21, -87636.05),
                           y = c(483099.5, 2194909.7)) #y = c(533099.5, 1894909.7))

#Alaska layer
ebs_layers <- akgfmaps::get_base_layers(select.region = "ebs", set.crs = "EPSG:3338")
ak_sppoly<-as(ebs_layers$akland, 'Spatial')

#create directory
dir.create('./shapefiles/',showWarnings = FALSE)

#name shapefiles 
shfiles<-c('EBSshelfThorson','NBSThorson','EBSslopeThorson')

#get files from google drive and set up
files<-googledrive::drive_find()
1 #for dvilasg@uw.edu

#get id shared folder from google drive
id.bering.folder<-files[which(files$name=='Shapefiles'),'id']

#list of files and folder
files.1<-googledrive::drive_ls(id.bering.folder$id)

#loop over shapefiles
for (i in shfiles) {
  
  #i=shfiles[1]
  
  id.data<-files.1[which(grepl(i,files.1$name)),]
  
  for (j in 1:nrow(id.data)) {
    
    #download data
    googledrive::drive_download(file=id.data$id[j],
                                path = paste0('./shapefiles/',id.data$name[j]),
                                overwrite = TRUE)
    
  }
  
  #shapefile EBS
  sh<-rgdal::readOGR(dsn='./shapefiles/',layer = i)
  
  if (i=='EBSslopeThorson') {
    
    #reproject shapefile
    proj4string(sh) <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0") 
    sh<-spTransform(sh,CRSobj = CRS('+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs'))
    
  }
  
  #shapefile name
  shname<-paste0(gsub('Thorson','',i),'_sh')
  
  #assign shapefiles
  assign(shname,sh)
  
}

#merge shapefiles
bs_sh1<-terra::union(EBSshelf_sh,NBS_sh)
bs_sh<-terra::union(bs_sh1,EBSslope_sh)
nbs_sh<-NBS_sh
ebs_sh<-EBSshelf_sh

#color palette
pal<-wesanderson::wes_palette('Zissou1',21,type='continuous')

##############################
#FIT PROJECT SETTINGS
#############################
#load fit file
load('./shelf EBS NBS VAST/Gadus macrocephalus/temp3d/b2_19822022fit.RData')

#add covariate data for the projected years into the fit$covariate_data
fit$covariate_data

#yrs
yrs<-as.integer(fit$year_labels)

#how manyt projected years we want
n_proj<-1

#project_yrs
project_yrs<-(last(yrs)+1):(last(yrs)+n_proj)

##############################
# 
#############################

df_proj<-subset(fit$covariate_data,Year==2022)
df_proj$BotTemp<-df_proj$BotTemp+0.5
df_proj$Year<- project_yrs

fit$covariate_data<-rbind(fit$covariate_data,df_proj)

#lets add BotTemp for the year 2011 
############################
#COVARIATE_DATA FROM DATA_GEOSTAT FILE
############################

#read data_geostat_temp file
df1<-readRDS(paste0('./data processed/Gadus macrocephalus/data_geostat_temp.rds'))
#df1[which(df1$year==2020),'bottom_temp_c']<-NA
df2<-subset(df1,year %in% project_yrs)

#select rows and rename
df3<-df2[,c("lat_start","lon_start","year",'scientific_name','weight_kg','effort','depth_m','LogDepth',"ScaleLogDepth",'Scalebottom_temp_c','bottom_temp_c','survey_name')]
colnames(df3)<-c('Lat','Lon','Year','Species','CPUE_kg','Effort','Depth','LogDepth','ScaleLogDepth','ScaleBotTemp','BotTemp','Region')

#data geostat
df4<-subset(df3,Region %in% c("Eastern Bering Sea Crab/Groundfish Bottom Trawl Survey",
                              "Northern Bering Sea Crab/Groundfish Survey - Eastern Bering Sea Shelf Survey Extension"))

#covariate data - filter by year and complete cases for env variables
#covariate_data<-subset(df2,Year>=yrs_region[1] & Year<=yrs_region[2])
new_covariate_data<-df3[complete.cases(df3[,c('BotTemp')]),] #,'ScaleLogDepth'
new_covariate_data1<-new_covariate_data[,names(fit$covariate_data)]
new_covariate_data1$CPUE_kg<-NA

#rbind with new year
fit$covariate_data<-rbind(fit$covariate_data,new_covariate_data1)

#project model example
p1<-project_model(x = fit,n_proj = n_proj)

summary(p1)
D_gt<-p1$D_gct[,1,]
D_gt[,paste0(project_yrs)]

############################
#COVARIATE_DATA FROM GRID WITH TEMP FROM ROMS
############################

#load fit file
load('./shelf EBS NBS VAST/Gadus macrocephalus/temp3d/b1_19822010fit.RData')

#read grid file with temp from ROMS
grid<-readRDS('./data processed/grid_slope_shelf_EBS_NBS_covariate_data.rds')
grid1<-grid[,c("Year","Lon","Lat","Temp","STRATA")]
grid2<-subset(grid1,STRATA!='EBSslope' & Year %in% project_yrs)

#create new df covariate data for projected years
new_covariate_data<-data.frame('Year'=grid2$Year,
                               'Lat'=grid2$Lat,
                               'Lon'=grid2$Lon,
                               "ScaleLogDepth"=NA,
                               "LogDepth"=NA,
                               "ScaleBotTemp"=NA,
                               "BotTemp"=grid2$Temp,
                               "CPUE_kg"=NA )

#rbind with new year
fit$covariate_data<-rbind(fit$covariate_data,new_covariate_data)

#project model example
p2<-project_model(x = fit,n_proj = n_proj)

############################
#COVARIATE_DATA FROM GRID WITH TEMP FROM ROMS but increasing temp
############################

#temperature increasing
add_temp<-1

#load fit file
load('./shelf EBS NBS VAST/Gadus macrocephalus/temp3d/b1_19822010fit.RData')

#read grid file with temp from ROMS
grid<-readRDS('./data processed/grid_slope_shelf_EBS_NBS_covariate_data.rds')
grid1<-grid[,c("Year","Lon","Lat","Temp","STRATA")]
grid2<-subset(grid1,STRATA!='EBSslope' & Year %in% project_yrs)

#create new df covariate data for projected years
new_covariate_data<-data.frame('Year'=grid2$Year,
                               'Lat'=grid2$Lat,
                               'Lon'=grid2$Lon,
                               "ScaleLogDepth"=NA,
                               "LogDepth"=NA,
                               "ScaleBotTemp"=NA,
                               "BotTemp"=grid2$Temp+add_temp,
                               "CPUE_kg"=NA )

#rbind with new year
fit$covariate_data<-rbind(fit$covariate_data,new_covariate_data)

#project model example
p3<-project_model(x = fit,n_proj = n_proj)


#loop
for (proj in c('p1','p2','p3')) {
  
  proj<-'p3'
  
  #get object
  p<-get(proj)
  
  #get predictions
  D_gt<-p$D_gct[,1,]
  #D_gt_proj<-D_gt[,paste0(project_yrs)]
  D_gt<-drop_units(D_gt)
  D_gt<-data.frame('cell'=c(1:fit$spatial_list$n_g),D_gt)
  D_gt<-D_gt[!duplicated(as.list(D_gt))]
  colnames(D_gt)<-c('cell',c(fit$year_labels,project_yrs))
  D_gt1<-reshape2::melt(D_gt,id=c('cell'))
  
  mdl <- make_map_info(Region = fit$settings$Region,
                       spatial_list = fit$spatial_list,
                       Extrapolation_List = fit$extrapolation_list)
  
  D <- merge(D_gt1, mdl$PlotDF, by.x='cell', by.y='x2i')
  
  # #extract with extrapolation mesh
  # mesh<-fit$spatial_list
  # class(fit$spatial_list)
  # plot(mesh$MeshList)
  # mesh<-fit$spatial_list
  # mesh1<-mesh$MeshList$anisotropic_mesh
  # class(mesh1)
  # 
  # #looking at the spatial field and what it looks like
  # gproj <- inla.mesh.projector(mesh,  dims = c(300, 300))
  # g.mean <- inla.mesh.project(gproj, model1$summary.random$spatial.field$mean)
  # g.sd <- inla.mesh.project(gproj, model1$summary.random$spatial.field$sd)
  # 
  #plot list to store plots
  plot_list<-list()
  
  #loop over years
  for (year in project_yrs) {
    
    year<-project_yrs[1]
    
    #subset by year
    D1<-subset(D,variable==year)
    D2<-D1[,c("value","Lat","Lon")]
    
    #df to spatialpoint df
    coordinates(D2) <- ~ Lon + Lat
    crs(D2)<-c(crs="+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
    
    #reproject coordinates for plotting purposes
    D2_1<-spTransform(D2,'+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs')
    D2_2<-data.frame(D2_1)
    
    #plot density map
    p<-ggplot() +
      geom_point(data=D2_2, aes(Lon, Lat, color=log(as.vector(value)), group=NULL),
                 ## These settings are necessary to avoid
                 ## overlplotting which is a problem here. May need
                 ## to be tweaked further.
                 size=1.2, stroke=0,shape=16)+ #tune size to remove blanks on the map, depending on the resolution of the tiff
      geom_polygon(data=ak_sppoly,aes(x=long,y=lat,group=group),fill = 'grey60')+
      scale_x_continuous(expand = c(0,0),
                         breaks = c(-175,-170,-165,-160))+
      scale_y_continuous(expand = c(0,0),
                         breaks = c(62,58,54))+
      geom_polygon(data=nbs_sh,aes(x=long,y=lat,group=group),fill=NA,col='black')+
      geom_polygon(data=ebs_sh,aes(x=long,y=lat,group=group),fill=NA,col='black')+
      coord_sf(crs = '+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs',
               xlim = panel_extent$x,
               ylim = panel_extent$y)+
      scale_color_gradientn(colours = pal,name=('log(kg/km²)'),
                            guide = guide_colorbar(  frame.colour = "black",ticks.colour = 'black'))+
      #annotation_north_arrow(location = "tr", which_north = "true",pad_x = unit(0.01, 'points'), pad_y = unit(10, 'points'),
      #                       style = north_arrow_fancy_orienteering(line_width = 1.5, text_size =6))+
      theme(panel.border = element_rect(fill=NA),aspect.ratio = 1,panel.grid.major = element_line(color = rgb(0, 0, 0,20, maxColorValue = 285), linetype = 'dashed', linewidth =  0.5),
            panel.background = element_rect(fill = NA),panel.ontop = TRUE,text = element_text(size=11),
            legend.background =  element_rect(fill = "transparent", colour = "transparent"),legend.key.height= unit(12, 'points'),
            legend.key.width= unit(12,units = 'points'),axis.title = element_blank(),legend.position = c(0.12,0.18),
            plot.title = element_text(size=16,vjust = -10, hjust=0.95),legend.title = element_text(size = 8),
            plot.margin = unit(c(-0.01,-0.01,-0.01,-0.01), "points"),legend.key = element_rect(color="black"))+
      #annotate("text", x = -256559, y = 1354909, label = "Alaska",parse=TRUE,size=7)+
      #annotate("text", x = -1176559, y = 1904909, label = "Russia",parse=TRUE,size=7)+
      labs(title=paste0(year))
    
    plot_list[[year]]<-p

}





