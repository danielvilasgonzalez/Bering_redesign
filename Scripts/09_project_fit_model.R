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
pack_cran<-c("splines",'ragg','ggplot2','cowplot')

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

#############################
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

#scenarios file
df_scn<-read.csv('./tables/SBT_scenarios.csv')
df_scn<-df_scn[,1:8]

#fit file
ff<-'fit.RData'

#create sp folder
dir.create('./output/species/')

#loop over species
for (sp in spp) {

  sp<-'Gadus macrocephalus'
  
  #create sp folder
  dir.create(paste0('./output/species/',sp))
  
##############################
# FIT PROJECT SETTINGS 
#############################

  #load fit file
  load(paste0('./shelf EBS NBS VAST/',sp,'/',ff))
  
  #add covariate data for the projected years into the fit$covariate_data
  #fit$covariate_data
  
  #yrs
  yrs<-as.integer(fit$year_labels)
  
  #how manyt projected years we want
  n_proj<-5
  
  #project_yrs
  project_yrs<-(last(yrs)+1):(last(yrs)+n_proj)

##############################
# PROJECT THROUGH SBT SCENARIOS
#############################

  #get scenarios
  df_scn<-read.csv('./tables/SBT_scenarios.csv')
  df_scn<-df_scn[,c(1:8)]
  
  #get raster stack
  stack_files<-list.files('./data processed/SBT scenarios/')
  
  #list covariate data for each scenario
  cov_list<-list()
  
  #list projected data for each scenario
  pr_list<-list()

    #loop over scenarios
    for (scn in unique(df_scn$scn_n)) {
      
      #scn<-unique(df_scn$scn_n)[2]
      
      #print scenario to check progress
      cat(paste(" #############     PROJECTING    #############\n",
                " #############   Species", sp, match(sp,spp), 'out of',length(spp),  "  #############\n",
                " #############  Scenario", scn, " #############\n"))

      #open stack of rasters
      st<-stack_files[grepl(paste0('scn',scn),stack_files)][1]
      st<-stack(paste0('./data processed/SBT scenarios/',st))
      
      #reproject shapefile
      proj4string(st) <- CRS('+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs') 
      
      #raster to points
      points<-data.frame(rasterToPoints(st))
      
      #create a df to store
      points3<-data.frame(matrix(nrow = 0,ncol = ncol(fit$covariate_data)))
      names(points3)<-names(fit$covariate_data)
      
      for (y in project_yrs) {
        
        #y<-project_yrs[1]
        
        #get points for year
        points1<-points[,c('x','y',paste0('y',y))]
        names(points1)<-c('Lon',"Lat",'BotTemp')
        
        #reproject df
        coordinates(points1)<- ~ Lon + Lat
        proj4string(points1) <- CRS('+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs') 
        points1<-spTransform(points1,CRSobj = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
        
        #create a new df
        points2<-data.frame(cbind(Year=y,Lat=points1$Lat,Lon=points1$Lon,ScaleLogDepth=NA,LogDepth=NA,ScaleBotTemp=NA,BotTemp=points1$BotTemp,CPUE_kg=NA))
        
        #add year
        points3<-rbind(points3,points2)
        
      }
      
      
      #add year to covariate data from fit
      cov_list[[scn]]<-points3
      
      #add to covariate data
      fit$covariate_data<-rbind(fit$covariate_data,points3)
    
      #project model example
      pm<-project_model(x = fit,n_proj = n_proj,)
      
      #add year to covariate data from fit
      pr_list[[paste0('scn',scn)]]<-pm
    }
  
  #save projection list
  save(pr_list, file = paste0("./output/species/",sp,'/fit_projection.RData'))
}

##############################
# PLOT PROJECTIONS
#############################

#create sp folder
dir.create('./figures/species/')

#loop over species
for (sp in spp) {
  
  sp<-'Gadus macrocephalus'
  
  #print scenario to check progress
  cat(paste(" #############      PLOTTING     #############\n",
            " #############   Species", sp, match(sp,spp), 'out of',length(spp),  "  #############\n",
            " #############  Scenario", scn, " #############\n"))
  
  #create sp folder
  dir.create(paste0('./figures/species/',sp))
  
  #save projection list
  load(paste0("./output/species/",sp,'/fit_projection.RData')) #pr_list
  
  for (proj in names(pr_list)) {
    
    #proj<-names(pr_list)[1]
    
    #get object
    p<-pr_list[[proj]]
    
    #load fit file
    load(paste0('./shelf EBS NBS VAST/',sp,'/',ff))
    
    #get predictions
    D_gt<-p$D_gct[,1,]
    #D_gt_proj<-D_gt[,paste0(project_yrs)]
    D_gt<-drop_units(D_gt)
    D_gt<-data.frame('cell'=c(1:fit$spatial_list$n_g),D_gt)
    
    colnames(D_gt)<-c('cell',c(fit$year_labels,project_yrs))
    D_gt1<-reshape2::melt(D_gt,id=c('cell'))
    
    mdl <- make_map_info(Region = fit$settings$Region,
                         spatial_list = fit$spatial_list,
                         Extrapolation_List = fit$extrapolation_list)
    
    D <- merge(D_gt1, mdl$PlotDF, by.x='cell', by.y='x2i')
    
    #plot list to store plots
    plot_list<-list()
    
    #loop over years
    for (y in project_yrs) {
      
      #y<-project_yrs[1]
      
      #subset by year
      D1<-subset(D,variable==y)
      D2<-D1[,c("value","Lat","Lon")]
      
      #df to spatialpoint df
      coordinates(D2) <- ~ Lon + Lat
      crs(D2)<-c(crs="+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
      
      #reproject coordinates for plotting purposes
      D2_1<-spTransform(D2,'+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs')
      D2_2<-data.frame(D2_1)
      
      #x and y cells
      xycells<-as.integer(sqrt(dim(D2_1)[1]))
      
      # create a template raster
      r1 <- raster(ext=extent(D2_1),ncol=xycells, nrow=xycells) #c(15800,15800) 7000
      
      #create raster
      r2<-rasterize(D2_1, r1 ,field='value')
      crs(r2) <- CRS('+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs')
      #plot(r2)
     
      #plot density map
      p<-
        ggplot() +
        geom_raster(data=as.data.frame(r2, xy = TRUE),aes(x=x,y=y,fill=log(layer)))+
        #geom_point(data=D2_2, aes(Lon, Lat, color=log(as.vector(value)), group=NULL),
                   ## These settings are necessary to avoid
                   ## overlplotting which is a problem here. May need
                   ## to be tweaked further.
                   #ize=1.2, stroke=0,shape=16)+ #tune size to remove blanks on the map, depending on the resolution of the tiff
        geom_polygon(data=ak_sppoly,aes(x=long,y=lat,group=group),fill = 'grey60')+
        scale_x_continuous(expand = c(0,0),
                           breaks = c(-175,-170,-165,-160))+
        scale_y_continuous(expand = c(0,0),
                           breaks = c(65,60,55))+
        geom_polygon(data=nbs_sh,aes(x=long,y=lat,group=group),fill=NA,col='black')+
        geom_polygon(data=ebs_sh,aes(x=long,y=lat,group=group),fill=NA,col='black')+
        coord_sf(crs = '+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs',
                 xlim = c(-1486559.21, -77636.05),
                 ylim = c(453099.5,2004909.7),
                 label_axes = "-NE-")+
        scale_fill_gradientn(colours = pal,name=('log(kg/km²)'),na.value = 'white',
                              guide = guide_colorbar(  frame.colour = "black",ticks.colour = 'black'))+
        #annotation_north_arrow(location = "tr", which_north = "true",pad_x = unit(0.01, 'points'), pad_y = unit(10, 'points'),
        #                       style = north_arrow_fancy_orienteering(line_width = 1.5, text_size =6))+
        theme(aspect.ratio = 1,panel.grid.major = element_line(color = rgb(0, 0, 0,20, maxColorValue = 285), linetype = 'dashed', linewidth =  0.5),
              panel.background = element_rect(fill = NA),panel.ontop = TRUE,text = element_text(size=10),
              legend.background =  element_rect(fill = "transparent", colour = "transparent"),legend.key.height= unit(10, 'points'),
              legend.key.width= unit(10, 'points'),axis.title = element_blank(),legend.position = c(0.12,0.23),
              panel.border = element_rect(fill = NA, colour = 'black'),legend.key = element_rect(color="black"),
              axis.text = element_text(color='black'),legend.spacing.y = unit(8, 'points'),
              axis.text.y.right = element_text(hjust= 0.1 ,margin = margin(-7,0,0,-25, unit = 'points'),color='black'),
              axis.text.x = element_text(vjust = 6, margin = margin(-7,0,0,-25, unit = 'points'),color='black'),
              axis.ticks.length = unit(-5,"points"),plot.title = element_text(size=12,vjust = -10, hjust=0.95,face="bold"),
              plot.margin=margin(c(5,5,5,5)),legend.title = element_text(size = 8))+
        #annotate("text", x = -256559, y = 1354909, label = "Alaska",parse=TRUE,size=7)+
        #annotate("text", x = -1176559, y = 1904909, label = "Russia",parse=TRUE,size=7)+
        labs(title=paste0(y))
      
      plot_list[[as.character(y)]]<-p
      
    }
    
    #scenario
    scn<-gsub('scn','',proj)
    
    #scenario name
    scn_name<-df_scn[which(df_scn$scn_n==scn),'Scenario']
    
    # now add the title
    title <- ggdraw() + 
      draw_label(
        paste0("Scn",scn,' - ',scn_name),
        fontface = 'bold',hjust=0.5,size=16,vjust=0.7
      ) 
    
    #save multiplot
    mp<-cowplot::plot_grid(plotlist = plot_list,nrow = 1,ncol = n_proj)
    ragg::agg_png(paste0('./figures/species/sp/projection_',proj,".png"), width = 20, height = 4.6, units = "in", res = 300)
    print(
        plot_grid(
        title, mp,
        ncol = 1,
        # rel_heights values control vertical title margins
        rel_heights = c(0.08, 1)),align = "hv")
    dev.off()
    
  }  
}  
