####################################################################
####################################################################
##
##    Project OM for each sp under multiple SBT conditions
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
version<-'VAST_v14_0_1'

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
2 #for dvilasg@uw.edu

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
load('./tables/SBT_projection.RData') #df_sbt

#load grid of NBS and EBS
load('./extrapolation grids/northern_bering_sea_grid.rda')
load('./extrapolation grids/eastern_bering_sea_grid.rda')
grid<-as.data.frame(rbind(data.frame(northern_bering_sea_grid,region='NBS'),data.frame(eastern_bering_sea_grid,region='EBS')))
grid$cell<-1:nrow(grid)

#number of simulation on model projection
n_sim<-100

#create sp folder
dir.create('./output/species/')

#loop over species
for (sp in spp) {

  #sp<-spp[1]
  
  #create sp folder
  dir.create(paste0('./output/species/',sp))
  
  #get list of fit data
  ff<-list.files(paste0('./shelf EBS NBS VAST/',sp),'fit',recursive = TRUE)
  
##############################
# FIT PROJECT SETTINGS 
#############################

  #load fit file
  load(paste0('./shelf EBS NBS VAST/',sp,'/',ff))
  #fit<-reload_model(fit)
  
  #check fit
  #check_fit(fit$parameter_estimates)
  #add covariate data for the projected years into the fit$covariate_data
  #fit$covariate_data
  
  #yrs
  yrs<-1982:2022
  
  #how manyt projected years we want
  n_proj<-5
  
  #project_yrs
  project_yrs<-(last(yrs)+1):(last(yrs)+n_proj)

##############################
# PROJECT THROUGH SBT SCENARIOS
#############################

  #get raster stack
  stack_files<-list.files('./data processed/SBT projections/')
  
  #list covariate data for each scenario
  #cov_list<-list()
  
  #list projected data for each scenario
  #pr_list<-list()

  #create folder simulation data
  dir.create(paste0('./output/species/',sp,'/simulated projected data/'))
  
    #loop over scenarios
    for (sbt in unique(df_sbt$sbt_n)) {
      
      #sbt<-unique(df_sbt$sbt_n)[1]
      
      #print scenario to check progress
      cat(paste(" #############     PROJECTING    #############\n",
                " #############   Species", sp, match(sp,spp), 'out of',length(spp),  "  #############\n",
                " #############  SBT", sbt, " #############\n"))

      #open stack of rasters
      st<-stack_files[grepl(paste0('SBT_',sbt),stack_files)][1]
      st<-stack(paste0('./data processed/SBT projections/',st))
      #plot(st)
      #title(main='x')
      
      #reproject shapefile
      #proj4string(st) <- CRS('+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs') 
      
      #raster to points
      points<-data.frame(rasterToPoints(st))
      
      #load fit file
      load(paste0('./shelf EBS NBS VAST/',sp,'/',ff))
      
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
      #cov_list[[sbt]]<-points3
      
      #add covariate data
      new_data<-rbind(fit$covariate_data,points3)
      
      #project model example
      pm<-VAST::project_model(x = fit,
                        n_proj = n_proj,
                        n_samples = n_sim,
                        new_covariate_data = new_data,
                        historical_uncertainty = 'none')
                        
      #historical_uncertainty = 'none')
    
      #ps<-pm
      
      #remove fit
      #rm(fit)
      #dyn.unload('C:/Program Files/R/R-4.2.2/library/VAST/executables/VAST_v13_1_0_TMBad.dll')
      
      #save list projections
      save(pm, file = paste0("./output/species/",sp,'/simulated projected data/fit_projection_SBT',sbt,'.RData'))
      
      rm(pm)
      gc()
      #add year to covariate data from fit
      #pr_list[[paste0('SBT',sbt)]]<-pm
    }
  
  #save projection list
  #save(pr_list, file = paste0("./output/species/",sp,'/simulated historical data/fit_projections.RData'))
}

##########################
# CHECK SIMULATED PROJECTIONS
##########################

#store results
df<-data.frame(matrix(NA,nrow = 0,ncol=4))
names(df)<-c('year','index','sim','sbt')

for (sbt in paste0('SBT',1:12)) {
  
  #print scenario to check progress
  cat(paste(" #############  ", sbt, " #############\n"))
  
  #load
  load(paste0('./output/species/Gadus macrocephalus/simulated projected data/fit_projection_',sbt,'.RData'))
  
  #check samples index fit
  #length(pm)
  
  #loop over 100 simulations
  for (sim in 1:length(pm)) {
    
    #sim<-1
    
    index<-
      df1<-data.frame(year=1982:2027,index=drop_units(pm[[sim]]$Index_ctl[,1:46,1])/1000,sim=sim,sbt=sbt)
    df<-rbind(df,df1)
  }
  
  rm(pm)
}

#mean projection
df2<-aggregate(index ~ year + sbt,df,FUN = function(x) c(mean = mean(x), sd = sd(x) ) )



###################################
# SBT projections
###################################

#save SBT table
load('./tables/SBT_projection.RData')#df_sbt

#name scenario
df_sbt$sbt<-paste0('SBT',df_sbt$sbt_n)
df_sbt$sbt2<-paste0(df_sbt$sbt,'_',df_sbt$Scenario)
df_sbt<-df_sbt[,c('sbt','sbt2')]

df1<-merge(df,df_sbt,by='sbt',all.x=TRUE)
df22<-merge(df2,df_sbt,by='sbt',all.x=TRUE)

df22$sbt2<-factor(df22$sbt2,levels = unique(df22$sbt2)[c(1,5:12,2:4)])
df1$sbt2<-factor(df1$sbt2,levels = unique(df22$sbt2)[c(1,5:12,2:4)])

#plot
ggplot()+
  geom_line(data=df1,aes(x=year,y=index,group=sim),color='grey80')+
  geom_ribbon(data=df22,aes(x=year,ymax=index[,'mean']+index[,'sd'],ymin=index[,'mean']-index[,'sd']),color='grey30')+
  geom_line(data=df22,aes(x=year,y=index[,'mean']),color='black')+
  theme_bw()+
  facet_wrap(~sbt2)

##############################
# PLOT MAP PROJECTIONS
#############################

#create sp folder
dir.create('./figures/species/')

#loop over species
for (sp in spp) {
  
  sp<-'Gadus macrocephalus'

  #create sp folder
  dir.create(paste0('./figures/species/',sp))
  
  #save projection list
  load(paste0("./output/species/",sp,'/fit_projection.RData')) #pr_list
  
  #load fit file
  load(paste0('./shelf EBS NBS VAST/',sp,'/',ff))
  
  # index_all<-array(NA,
  #                  dim=list(length(names(pr_list)),length(1982:2027),length(spp)),
  #                  dimnames = list(names(pr_list),c(1982:2027),spp))
  
  for (proj in names(pr_list)) {
    
    #proj<-names(pr_list)[1]
    
    #print scenario to check progress
    cat(paste(" #############      PLOTTING     #############\n",
              " #############   Species", sp, match(sp,spp), 'out of',length(spp),  "  #############\n",
              " #############  ", proj, " #############\n"))
    
    
    #get object
    p<-pr_list[[proj]]
    
  #   sum(p$Index_gctl)
  #   
  #   index<-p$Index_gctl[,1,,]
  #   i<-colSums(index)
  #   
  #   index_all[proj,,sp]<-i
  # }
  #   
  # x<-reshape2::melt(index_all[,,sp])
  # 
  # ggplot()+
  #   geom_line(data=x,aes(x=Var2,y=value,color=Var1))+
  #   facet_wrap(~Var1)

    
    #get predictions
    D_gt<-p$D_gct[,1,]
    #D_gt_proj<-D_gt[,paste0(project_yrs)]
    D_gt<-drop_units(D_gt)
    D_gt<-data.frame('cell'=c(1:p$n_g),D_gt)
    
    #p$
    
    colnames(D_gt)<-c('cell',c(fit$year_labels,project_yrs))
    D_gt1<-reshape2::melt(D_gt,id=c('cell'))
    
    mdl <- make_map_info(Region = fit$settings$Region,
                         spatial_list = fit$spatial_list,
                         Extrapolation_List = fit$extrapolation_list)
    
    D <- merge(D_gt1, mdl$PlotDF, by.x='cell', by.y='x2i')
    
    #plot list to store plots
    plot_list<-list()
    
    #create stack
    proj_stack<-stack()
    
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
     
      proj_stack<-raster::addLayer(proj_stack,r2)
      
      #plot density map
      p<-
        ggplot() +
        geom_raster(data=as.data.frame(r2, xy = TRUE),aes(x=x,y=y,fill=log(layer)))+
        #geom_point(data=D2_2, aes(Lon, Lat, color=log(as.vector(value)), group=NULL),
                   ## These settings are necessary to avoid
                   ## overlplotting which is a problem here. May need
                   ## to be tweaked further.
                   #size=1.2, stroke=0,shape=16)+ #tune size to remove blanks on the map, depending on the resolution of the tiff
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
    
    #save raster stack
    writeRaster(proj_stack, paste0('output/species/',sp,'/simulated projected data/biomass_projection_',proj,'.grd'),overwrite=TRUE)
    
    #scenario
    #sbt<-gsub('sbt','SBT',proj)
    
    #scenario name
    #sbt_name<-df_sbt[which(grepl(paste0(sbt,'\\>'),paste0('SBT',df_sbt$sbt_n))),'Scenario']
    sbt_name<-df_sbt[which(grepl(paste0(proj,'\\>'),df_sbt$sbt)),'Scenario']
    
    # now add the title
    title <- ggdraw() + 
      draw_label(
        paste0(proj,' - ',sbt_name),
        fontface = 'bold',hjust=0.5,size=16,vjust=0.7
      ) 
    
    #save multiplot
    mp<-cowplot::plot_grid(plotlist = plot_list,nrow = 1,ncol = n_proj)
    ragg::agg_png(paste0('./figures/species/',sp,'/projected_distribution_',proj,".png"), width = 20, height = 4.6, units = "in", res = 300)
    print(
        plot_grid(
        title, mp,
        ncol = 1,
        # rel_heights values control vertical title margins
        rel_heights = c(0.08, 1)),align = "hv")
    dev.off()
    
  }  
}  

###################################
# CHECK INDICES
##################################

#load grid of NBS and EBS
load('./extrapolation grids/northern_bering_sea_grid.rda')
load('./extrapolation grids/eastern_bering_sea_grid.rda')
grid<-as.data.frame(rbind(data.frame(northern_bering_sea_grid,region='NBS'),data.frame(eastern_bering_sea_grid,region='EBS')))
grid$cell<-1:nrow(grid)

#true
load(paste0('./shelf EBS NBS VAST/',sp,'/',ff)) #fit
true_ind<-data.frame('year'= as.integer(names(fit$Report$Index_ctl[,,1])),
                     'value'=as.vector(fit$Report$Index_ctl[,,1]),
                     'sbt'='true')



dens<-fit$Report$D_gct[,1,]
area<-grid$Area_in_survey_km2
bio<-sweep(dens, 2, area, FUN="*")
bio_df<-data.frame(bio_t=colSums(bio)/1000,year=1982:2022)



#projection list
 load( file = paste0("./output/species/",sp,'/fit_projection.RData')) #pr_list
 ind<-data.frame(matrix(nrow = 0,ncol=3))
 colnames(ind)<-c('year','value','sbt')
 ind_area<-ind
 
 for (i in names(pr_list)) {
   
   #i<-names(pr_list)[1]
   pr<-pr_list[[i]]
   df<-pr['Index_ctl']
   df1<-data.frame(year=1982:2027,
                   value=as.vector(df$Index_ctl[1,,1]),
                   sbt=i)
   ind<-rbind(ind,df1)
 
   dens_pr<-pr$D_gct[,1,]
   bio_pr<-sweep(dens_pr, 2, area, FUN="*")
   df1_pr<-data.frame(year=1982:2027,
                      value=as.vector(colSums(bio_pr)),
                      sbt=i)
   
   ind_area<-rbind(ind_area,df1_pr)
 }

 ind_area$output<-'sum(dens*area)'
 ind$output<-'Index_fit'
 ind1<-rbind(ind_area,ind)
 df_sbt$sbt<-paste0('SBT',df_sbt$sbt_n)
   
 ind2<-merge(ind1,df_sbt,by='sbt',all.x=TRUE)
 dim(ind1)
 dim(ind2)
 
 ind2$match<-ifelse(ind2$sbt %in% c('SBT2','SBT8','SBT10'),'non','yes')
 ind2$sbt2<-paste(ind2$sbt,ind2$Scenario)
 ind2$sbt2<-factor(ind2$sbt2,levels = unique(ind2$sbt2)[c(1,5:12,2:4)])
 
 #plot comparison
 ggplot()+
   #geom_line(data=bio_df,aes(x=1982:2022,y=bio_t),color='black',linetype='dashed')+
   geom_line(data=subset(ind2,output=='Index_fit'),aes(x=year,y=value/1000,color=sbt2,linetype=sbt2),alpha=0.8)+
   #ggrepel::geom_text_repel(data = subset(ind1, year == "2022" & output =='Index_fit'), 
  #                          aes(label = sbt, colour = sbt, x = 2022, y = value/1000), 
  #                          box.padding = 0.5, max.overlaps = Inf) +
   #geom_line(data=ind_area,aes(x=year,y=value,color=sbt),linetype='dashed',alpha=0.8)+
   scale_linetype_manual(values=c('solid','dashed',rep('solid',5),'dashed','solid','dashed','solid','solid'),name='SBT projections')+
   scale_color_manual(values = hue_pal()(12),name='SBT projections')+
   geom_line(data=true_ind,aes(x=year,y=value/1000),color='black')+
   labs(y='t',x='year')+
   scale_x_continuous(breaks=c(1982:2027),limits = c(1982,2027))+
   theme_bw()+
   #guides(linetype = 'none')+
   theme(axis.text.x = element_text(angle=90,vjust = 0.5),panel.grid.minor.x = element_blank())
 
 # #plot comparison zoom in
 # ggplot()+
 #   geom_line(data=ind,aes(x=year,y=value/1000,color=sbt),linewidth=1,alpha=0.8)+
 #   #geom_point(data=true_ind,aes(x=year,y=value/1000,fill=sbt))+
 #   labs(y='t',color='projection',fill='fitted')+
 #   scale_x_continuous(breaks=c(1982:2027),limits = c(2023,2027))+
 #   scale_y_continuous(limits = c(300000,1000000))+
 #   scale_color_discrete(labels=paste0('sbt',1:9,'_',c("status quo","gradually cold","gradually warm","medium variation I","medium variation II","high variation I","high variation II",   
 #                               "extreme variation I","extreme variation II")))+
 #   theme_bw()+
 #   theme(axis.text.x = element_text(angle=90,vjust = 0.5),panel.grid.minor.x = element_blank())
 #  
