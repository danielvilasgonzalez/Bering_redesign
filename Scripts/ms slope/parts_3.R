#############################################
#
#
#   Evaluate if changes on the temp have caused changes on the bathymetric distribution of species
#   (warm years tend to move species north or deeper)
#
#
##############################################

#####################################
# Settings
#####################################

#clear all objects
rm(list = ls(all.names = TRUE)) 
#free up memrory and report the memory usage
gc() 

#libraries from cran to call or install/load
pack_cran<-c('cowplot','ggspatial','raster','rasterVis','rgeos','scales','rnaturalearth','grid','ggplot2','lubridate','ragg','rgdal','dplyr','sp','ggtext','magick','shadowtext')

#install pacman to use p_load function - call library and if not installed, then install
if (!('pacman' %in% installed.packages())) {
  install.packages("pacman")}

#install akgfmaps to extract shapefile of Alaska
if (!('akgfmaps' %in% installed.packages())) {
  devtools::install_github('afsc-gap-products/akgfmaps')};library(akgfmaps)

#install VAST
if (!('VAST' %in% installed.packages())) {
  devtools::install_github("james-thorson/VAST@main", INSTALL_opts="--no-staged-install")};library(VAST)

#load/install packages
pacman::p_load(pack_cran,character.only = TRUE)

#set working directory
#out_dir<-'C:/Users/Daniel.Vilas/Work//Adapting Monitoring to a Changing Seascape/'
out_dir<-'/Users/daniel/Work/Adapting Monitoring to a Changing Seascape/'
setwd(out_dir)

####################
# GRIDS
####################

#load grid of NBS and EBS
load('./extrapolation grids/northern_bering_sea_grid.rda')
load('./extrapolation grids/eastern_bering_sea_grid.rda')
load('./extrapolation grids/bering_sea_slope_grid.rda')

# #nbs + ebs grid
# nbs_grid<-data.frame(northern_bering_sea_grid[,c('Lat','Lon','Area_in_survey_km2')])
# nbs_grid$region<-'NBS'
# ebs_grid<-data.frame(eastern_bering_sea_grid[,c('Lat','Lon','Area_in_survey_km2')])
# ebs_grid$region<-'EBSshelf'
# nbsebs_grid<-rbind(nbs_grid,ebs_grid)
# nbsebs_grid$Site<-1:nrow(nbsebs_grid)
# 
# #slope grid
# slope_grid<-data.frame(bering_sea_slope_grid[,c('Lat','Lon','Area_in_survey_km2')])
# slope_grid$region<-'EBSslope'
# slope_grid$Site<-1:nrow(slope_grid)
# 
# #all grid
# all_grid<-rbind(nbsebs_grid,slope_grid)

####################
# convergence table and converged species for EBS+NBS and EBSslope models
####################

conv_tab<-read.csv('./tables/slope_ebsnbs_convspp.csv')
sp_shelfslope<-conv_tab[which(conv_tab$slope=='There is no evidence that the model is not converged' & conv_tab$EBS_NBS=='There is no evidence that the model is not converged'),'spp']
sp_shelf<-conv_tab[which(conv_tab$EBS_NBS=='There is no evidence that the model is not converged'),'spp']
sp_slope<-conv_tab[which(conv_tab$slope=='There is no evidence that the model is not converged' ),'spp']

#####################################
# Polygon land
#####################################

# Define plot extent (through trial end error) units km
panel_extent <- data.frame(x = c(-1716559.21, -77636.05), #x = c(-1326559.21, -87636.05),
                           y = c(483099.5, 2194909.7)) #y = c(533099.5, 1894909.7))


#Alaska land shapefile layer
ebs_layers <- akgfmaps::get_base_layers(select.region = "ebs", set.crs = "EPSG:3338")
ak_sppoly<-as(ebs_layers$akland, 'Spatial')

#####################################
# Polygon regions shapefiles (EBS, NBS and slope)
#####################################

#name shapefiles 
shfiles<-c('EBSshelfThorson','NBSThorson','EBSslopeThorson')

#loop over shapefiles
for (i in shfiles) {
  
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
bs_sh1<-raster::union(EBSshelf_sh,NBS_sh)
bs_sh<-raster::union(bs_sh1,EBSslope_sh)

####################
# average SBT by region
####################

#get grid
load(file = './data processed/grid_EBS_NBS.RData')

#get sites grid
grid.ebs_year_sites<-unique(grid.ebs_year[,c("DepthGEBCO","region","Area_in_survey_km2","Lat","Lon")])
grid.ebs_year_sites$Site<-c(1:nrow(subset(grid.ebs_year_sites,region %in% c('NBS', 'EBSshelf'))),1:nrow(subset(grid.ebs_year_sites,region %in% c('EBSslope'))))
grid.ebs_year_sitesslope<-grid.ebs_year_sites[which(grid.ebs_year_sites$region=='EBSslope'),]
grid.ebs_year_sitesslope$region[grid.ebs_year_sitesslope$Lat>57.5]<-'upperEBSslope'
grid.ebs_year_sitesslope$region[grid.ebs_year_sitesslope$Lat<57.5]<-'lowerEBSslope'

#get detph in grid
grid.ebs_year1<-unique(grid.ebs_year[,c("DepthGEBCO","region","Area_in_survey_km2","Lat","Lon")])
nrow(grid.ebs_year1[which(grid.ebs_year1$region=='NBS'),])
#grid.ebs_year1s<-unique(grid.ebs_year[,c('Stratum',"DepthGEBCO","region","Area_in_survey_km2","Lat","Lon")])

# #add depths regions
# region_depth<-ifelse(grid.ebs_year1$DepthGEBCO<100,'0_100',
#                      ifelse(grid.ebs_year1$DepthGEBCO>=100 & grid.ebs_year1$DepthGEBCO<=200,'100_200',
#                             ifelse(grid.ebs_year1$DepthGEBCO>=200 & grid.ebs_year1$DepthGEBCO<=300,'200_300',
#                                    ifelse(grid.ebs_year1$DepthGEBCO>=300 & grid.ebs_year1$DepthGEBCO<=400,'300_400',NA))))
# 
# region_depth<-ifelse(grid.ebs_year1$region=='NBS',NA,region_depth)

# #cbind
# grid.ebs_year1<-cbind(grid.ebs_year1,region_depth)

#plot
ggplot()+
  geom_point(data=grid.ebs_year1,aes(x=Lon,y=Lat,color=region))

# #average depth by region and year
# depth_region<-
#   aggregate(DepthGEBCO ~ region + region,grid.ebs_year,FUN=mean)
# grid.ebs_year_d1<-subset(grid.ebs_year,region %in% c('EBSshelf','NBS'))
# depth_region<-
#   rbind(depth_region,
#         c('EBS+NBS',mean(grid.ebs_year_d1$DepthGEBCO)),
#         c('all',mean(grid.ebs_year$DepthGEBCO)))


#average SBT by region and year
temp_region<-
  aggregate(Temp ~ Year + region,grid.ebs_year,FUN=mean)
names(temp_region)<-c('year','region','temp')

# ebsnbstemp<-
#   aggregate(Temp ~ Year ,subset(grid.ebs_year,region %in% c('EBSshelf','NBS')),FUN=mean)
# ebsnbstemp<-data.frame('year'=ebsnbstemp$Year,
#                        'region'='EBS+NBS',
#                        'temp'=ebsnbstemp$Temp)

#get average SBT for all region (EBSshelf+NBS+EBSslope)
alltemp<-aggregate(temp ~ year, temp_region,FUN=mean)
alltemp$region<-'all'
alltemp<-alltemp[,c("year","region","temp")]

#rbind
temp_region<-rbind(temp_region,alltemp)

#plot
#p<-
  ggplot()+
  #geom_hline(yintercept = mean(alltemp$temp),color='black')+
  geom_hline(yintercept = mean(subset(alltemp,year %in% c(2002:2016))$temp),color='grey30')+
  geom_rect(aes(xmin=2002,xmax=2016,ymin=-Inf,ymax=Inf),alpha=0.2)+
  geom_vline(xintercept = c(2003:2005,2015,2016),linetype='dashed',color='red')+
  geom_vline(xintercept = c(2008,2009,2012),linetype='dashed',color='blue')+
  geom_point(data=temp_region,aes(x=year,y=temp,color=region))+
  scale_color_manual(values=c('EBSslope'='#B4AF46','EBSshelf'='#B4464B','NBS'='#4682B4',all='black'),
                     breaks = c('all','EBS+NBS','EBSshelf','NBS','EBSslope'),
                     labels=c('ALL','EBS+NBS','EBSshelf','NBS','EBSslope'))+
  scale_x_continuous(breaks=c(seq(from=1982,to=2022,by=5)))+
  theme_bw()+
  labs(x='',y='average SBT')

#save  plot
agg_png(paste0('./figures slope/SBT_regions.png'), width = 7, height = 4, units = "in", res = 300)
print(p)
dev.off()

# Scaling values by group
temp_region1 <- data.frame(temp_region %>%
                             group_by(region) %>%
                             mutate(scaledtemp = scale(temp)) %>%
                             ungroup())

#classify by cold/warm/normal years
temp_region1$typeyear <- ifelse(temp_region1$scaledtemp > 1, "warm", ifelse(temp_region1$scaledtemp < -1, "cold", "normal"))

#remove slope
temp_region11<-subset(temp_region1,region %in% c('EBSshelf','all','NBS'))

#cold/warm years
cyrs<-unique(temp_region1[which(temp_region1$typeyear == c('cold') & temp_region1$region=='all'),'year'])
wyrs<-unique(temp_region1[which(temp_region1$typeyear == c('warm') & temp_region1$region=='all'),'year'])

# Define the years you want to bold
bold_years <- c(2002, 2004, 2008, 2010, 2012, 2016)

# Create custom labels function
custom_labels <- function(x) {
  sapply(x, function(year) {
    if (year %in% bold_years) {
      paste0("**", year, "**")
    } else {
      as.character(year)
    }
  })
}

#plot
p<-
ggplot() +
  geom_rect(aes(xmin = 1982, xmax = 2022, ymin = -1, ymax = 1), linetype = 'dashed', alpha = 0.2) +
  geom_rect(aes(xmin = 2002, xmax = 2016, ymin = -Inf, ymax = Inf), alpha = 0.2) +
  geom_vline(xintercept = cyrs, linetype = 'dashed', color = 'blue') +
  geom_vline(xintercept = wyrs, linetype = 'dashed', color = 'red') +
  geom_hline(yintercept = 0, alpha = 0.5, linetype = 'dashed') +
  geom_point(data = temp_region1, aes(x = year, y = scaledtemp, color = region)) +
  geom_line(data = subset(temp_region1, region == 'all'), aes(x = year, y = scaledtemp, color = region)) +
  scale_color_manual(values = c('EBSslope' = '#B4AF46', 'EBSshelf' = '#B4464B', 'NBS' = '#4682B4', 'EBS+NBS' = '#b649b1', 'all' = 'black'),
                     breaks = c('all', 'EBS+NBS', 'EBSshelf', 'NBS', 'EBSslope'),
                     labels = c('ALL', 'EBS+NBS', 'EBSshelf', 'NBS', 'EBSslope')) +
  theme_bw() +
  scale_x_continuous(breaks = c(1985,1990,1995,2000,bold_years,2020), labels = custom_labels, limits = c(1982, 2022)) +
  theme(axis.text.x = element_markdown(angle = 45,hjust=1),
        axis.text.y = element_text()) +
  labs(x = '', y = 'average scaled SBT')

#geom_vline(xintercept = c(2008,2009,2012),linetype='dashed',color='blue',alpha=0.5)

#save  plot
agg_png(paste0('./figures slope/scaledSBT_regions1.png'), width = 7, height = 4, units = "in", res = 300)
print(p)
dev.off()

#temp_region1[which(temp_region1$typeyear %in% c('warm','cold')),'year']

####################
# average SBT by region
####################

for (sp in sp_shelfslope) {
  
  sp<-sp_shelfslope[2]
  
  #sp<-'Gadus macrocephalus'
  
  #read data_geostat_temp file
  df1<-readRDS(paste0('./data processed/species/',sp,'/data_geostat_temp.rds'))

  #for slope data
  df11<-subset(df1,survey_name== "Eastern Bering Sea Slope Bottom Trawl Survey")
  #df11$survey_name[df11$survey_name == 'Eastern Bering Sea Slope Bottom Trawl Survey'] <- 'EBSslope'
  df11$survey_name[df11$survey_name == 'Eastern Bering Sea Slope Bottom Trawl Survey' & df11$lat_start >57.5] <- 'upperEBSslope'
  df11$survey_name[df11$survey_name == 'Eastern Bering Sea Slope Bottom Trawl Survey' & df11$lat_start <57.5] <- 'lowerEBSslope'
  
  #yrs only for slope
  yrs<-unique(df11$year)
  
  #for shelf data
  df12<-subset(df1,survey_name== "Eastern Bering Sea Crab/Groundfish Bottom Trawl Survey")
  df12$survey_name[df12$survey_name == 'Eastern Bering Sea Crab/Groundfish Bottom Trawl Survey'] <- 'EBSshelf'
  
  #for shelf data deeper than 106 (3rd quartile)
  df13<-subset(df1,survey_name== "Northern Bering Sea Crab/Groundfish Survey - Eastern Bering Sea Shelf Survey Extension")
  df13$survey_name[df13$survey_name == 'Northern Bering Sea Crab/Groundfish Survey - Eastern Bering Sea Shelf Survey Extension'] <- 'NBS'
  
  #rbind region specific df
  df1<-rbind(df11,df12,df13)
  
  # Define a data frame for the rectangles
  rect_data <- data.frame(
    xmin = c(2001.5, 2005.5, 2013.5),
    xmax = c(2005.5, 2013.5, 2022.5),
    fill = factor(c("warm", "cold", "warm"))
  )
  
  df1$survey_name <- factor(df1$survey_name, levels = c("NBS", "upperEBSslope","lowerEBSslope","EBSshelf", "warm", "cold"))
  
  
  # Create the plot
  p1<-
  ggplot() +
    geom_rect(data = rect_data, aes(ymin = 0, ymax = Inf, xmin = xmin, xmax = xmax, fill = fill), alpha = 0.3) +
    geom_boxplot(data = df1, aes(x = year, y = cpue_kgkm2, fill = survey_name, group = interaction(year, survey_name)), outliers = FALSE,position=position_dodge(width=0.8,preserve = "single")) +
    scale_x_continuous(limits = c(2001.5, 2022.5),breaks=c(2002:2022)) +
    theme_minimal()+
    theme(text = element_text(size=14),legend.position = 'none')+
    labs(title='OBSERVED DATA',x='')+
    scale_fill_manual(values = c("EBSshelf" = "green",'NBS'='purple', "lowerEBSslope" = "yellow1","upperEBSslope" = "yellow4", "warm" = "red", "cold" = "blue"), name = "env . region")
    #scale_fill_manual(values = c("EBSshelf" = "green",'NBS'='purple', "EBSslope" = "yellow", "warm" = "red", "cold" = "blue"), name = "env . region")
  
  #EBS+NBS model fit and density
  load(paste0('./shelf EBS NBS VAST/',sp,'/fit.RData'))
  fit_shelf<-fit
  dens_shelf<-drop_units(fit_shelf$Report$D_gct[,1,]) #dens
  dens_shelfnbs<-dens_shelf[1:15180,]
  dens_shelfebs<-dens_shelf[15181:nrow(dens_shelf),]  
  dens_shelf1nbs<-reshape2::melt(as.matrix(dens_shelfnbs))
  dens_shelf1ebs<-reshape2::melt(as.matrix(dens_shelfebs))
  dens_shelf1ebs$region<-'EBSshelf'
  dens_shelf1nbs$region<-'NBS'
  dens_shelf1<-rbind(dens_shelf1nbs,dens_shelf1ebs)

  #slope model fit ande density
  #load(paste0('./slope EBS VAST/',sp,'/fit.RData'))
  load(paste0('./slope EBS VAST/',sp,'/fit_st.RData'))
  fit_slope<-fit
  dens_slope<-drop_units(fit_slope$Report$D_gct[,1,]) #dens
  dens_slope<-cbind('Site'=1:3041,as.data.frame(dens_slope),'region'=grid.ebs_year_sitesslope$region)
  dens_slope1<-reshape2::melt(dens_slope,id.vars = c('Site','region'))
  dens_slope1<-dens_slope1[c('Site','variable','value','region')]
  names(dens_slope1)[2]<-'Time'
  dens_slope_catch<-dens_slope1
  dens_slope_catch$value<-dens_slope_catch$value*2
    
  #merge dens data
  dens_all<-rbind(dens_shelf1,dens_slope1)
  dens_all$Time<-as.numeric(dens_all$Time)
  dens_all$region <- factor(dens_all$region, levels = c("EBSshelf", "NBS", "upperEBSslope","lowerEBSslope", "warm", "cold"))
  
  #merge dens data with catch *2 for slope
  dens_all_catch<-rbind(dens_shelf1,dens_slope_catch)
  
  # Create the plot
  p2<-
  ggplot() +
    geom_rect(data = rect_data, aes(ymin = 0, ymax = Inf, xmin = xmin, xmax = xmax, fill = fill), alpha = 0.3) +
    scale_x_continuous(limits = c(2001.5, 2022.5), breaks = c(2002:2022)) +
    geom_boxplot(data = dens_all, aes(x = Time, y = value, fill = region, group = interaction(Time, region)), outliers = FALSE,position=position_dodge(width=0.8,preserve = "single")) +
    #geom_boxplot(data = subset(dens_all,region=='EBSslope'), aes(x = Time, y = value, fill = region, group = interaction(Time, region)), outliers = FALSE,position=position_dodge(width=0.8,preserve = "single")) +
    theme_minimal() +
    theme(text = element_text(size=14),legend.position = 'none')+
    labs(title='PREDICTED DATA',y='cpue_kgkm2',x='')+
    scale_fill_manual(values = c("EBSshelf" = "green",'NBS'='purple', "lowerEBSslope" = "yellow1","upperEBSslope" = "yellow4", "warm" = "red", "cold" = "blue"), name = "env . region")
  #scale_fill_manual(values = c("EBSshelf" = "green",'NBS'='purple', "EBSslope" = "yellow", "warm" = "red", "cold" = "blue"), name = "env . region")
  
  # Create the plot
  p3<-
    ggplot() +
    geom_rect(data = rect_data, aes(ymin = 0, ymax = Inf, xmin = xmin, xmax = xmax, fill = fill), alpha = 0.3) +
    scale_x_continuous(limits = c(2001.5, 2022.5), breaks = c(2002:2022)) +
    geom_boxplot(data = dens_all_catch, aes(x = Time, y = value, fill = region, group = interaction(Time, region)), outliers = FALSE,position=position_dodge(width=0.8,preserve = "single")) +
    #geom_boxplot(data = subset(dens_all,region=='EBSslope'), aes(x = Time, y = value, fill = region, group = interaction(Time, region)), outliers = FALSE,position=position_dodge(width=0.8,preserve = "single")) +
    theme_minimal() +
    theme(text = element_text(size=14),legend.position = 'none')+
    labs(title='PREDICTED DATA',y='cpue_kgkm2 (catch slope x2)',x='')+
    scale_fill_manual(values = c("EBSshelf" = "green",'NBS'='purple', "lowerEBSslope" = "yellow1","upperEBSslope" = "yellow4", "warm" = "red", "cold" = "blue"), name = "env . region")
  #scale_fill_manual(values = c("EBSshelf" = "green",'NBS'='purple', "EBSslope" = "yellow", "warm" = "red", "cold" = "blue"), name = "env . region")
  
  # Create the plot
  p4<-
    ggplot() +
    geom_rect(data = rect_data, aes(ymin = -Inf, ymax = Inf, xmin = xmin, xmax = xmax, fill = fill), alpha = 0.3) +
    scale_x_continuous(limits = c(2001.5, 2022.5), breaks = c(2002:2022)) +
    geom_boxplot(data = subset(dens_all,region %in% c('upperEBSslope','lowerEBSslope','NBS')), aes(x = Time, y = log(value), fill = region, group = interaction(Time, region)), outliers = FALSE, position=position_dodge(width=0.8,preserve = "single")) +
    #geom_boxplot(data = subset(dens_all,region=='EBSslope'), aes(x = Time, y = value, fill = region, group = interaction(Time, region)), outliers = FALSE,position=position_dodge(width=0.8,preserve = "single")) +
    theme_minimal() +
    theme(text = element_text(size=14),legend.position = 'none')+
    labs(title=paste0('SLOPE - PREDICTED DATA'),y='log(cpue_kgkm2)',x='')+
    scale_fill_manual(values = c("EBSshelf" = "green",'NBS'='purple', "lowerEBSslope" = "yellow1","upperEBSslope" = "yellow4", "warm" = "red", "cold" = "blue"), name = "env . region")
  #scale_fill_manual(values = c("EBSshelf" = "green",'NBS'='purple', "EBSslope" = "yellow", "warm" = "red", "cold" = "blue"), name = "env . region")
  
  #get legend
  legend1 <- cowplot::get_legend( 
    p1 + 
      theme(legend.position = "right",legend.spacing.y = unit(-0.1,'in')) #,legend.justification = c(0.72,0.5)
  ) 
  
  #plots
  plots<-plot_grid(p1,p2,p4,ncol=1)
  plots1<-
    plot_grid(plots,legend1,nrow=1,rel_widths = c(1,0.08))
  
  #title
  title <- ggdraw() + draw_label(sp, fontface='bold')
  
  #print plots
  print(
    plot_grid(title,plots1,ncol=1, rel_heights=c(0.05, 1))
  )
  
  ########################################################################
  ########################################################################
  
  # dens_all_catch1<-merge(dens_all_catch,grid.ebs_year_sites,by.x=c('Site','region'),by.y=c('Site','region'))
  # dens_all_catch1$bio<-dens_all_catch1$value*dens_all_catch1$Area_in_survey_km2
  # 
  # dens_all_catch2<-
  #   dens_all_catch1 %>% group_by(Time, region) %>%
  #   mutate(ind = sum(bio))
  # dens_all_catch3<-unique(dens_all_catch2[,c("Time","region","ind")])
  # 
  # dens_all_catch4<-
  #   dens_all_catch3 %>% group_by(Time) %>%
  #   mutate(frac = ind / sum(ind))
  # 
  # dens_all_catch4$Time<-as.factor(dens_all_catch4$Time)
  # 
  # # Your original plot code with correction for coloring axis text
  # #p<-
  #   ggplot() +
  #   geom_bar(data = dens_all_catch4, aes(x = Time, y = frac, fill = region), stat = 'identity') +
  #   scale_fill_manual(values = c('EBSslope' = '#B4AF46', 'EBSshelf' = '#B4464B', 'NBS' = '#4682B4'),
  #                     breaks = c('EBSshelf', 'NBS', 'EBSslope'),
  #                     labels = c('EBSshelf', 'NBS', 'EBSslope')) +
  #   #scale_x_discrete(labels = ind_df2$year) +
  #   labs(y='abundance proportion')+
  #   theme_bw() +
  #   theme(text=element_text(size=14))+
  #   #facet_wrap(~sp,scales='free_x', nrow = 5) +
  #   theme(axis.text.x = element_text(color = ifelse(levels(dens_all_catch4$Time) %in% c(2003:2005,2015:2016), "red", 
  #                                                   ifelse(levels(dens_all_catch4$Time) %in% c(2009,2008,2012), "blue",'black'))))
  # 
  
  ########################################################################
  ########################################################################
  
  ########################
  ########################
  # depth bins
  ########################
  ########################
  
  # Define a data frame for the rectangles
  rect_data <- data.frame(
    xmin = c(2001.5, 2005.5, 2013.5),
    xmax = c(2005.5, 2013.5, 2016.5),
    fill = factor(c("warm", "cold", "warm"))
  )
  
  dens_all1<-merge(dens_all,grid.ebs_year_sites,by.x=c('Site','region'),by.y=c('Site','region'))
  dens_all1<-dens_all1[which(dens_all1$DepthGEBCO>0),]
  dens_all1$depth_bin<-ifelse(dens_all1$DepthGEBCO<50,'0-50',
                         ifelse(dens_all1$DepthGEBCO<100,'50-100',
                                ifelse(dens_all1$DepthGEBCO<150,'100-150',
                                       ifelse(dens_all1$DepthGEBCO<200,'150-200','>200'))))
  
  dens_all1$depth_bin<-factor(dens_all1$depth_bin, levels = c("0-50", "50-100", '100-150',"150-200",">200"))

  ggplot()+
    geom_rect(data = rect_data, aes(ymin = -Inf, ymax = Inf, xmin = xmin, xmax = xmax, fill = fill), alpha = 0.3) +
    geom_boxplot(data=dens_all1,aes(x=Time,y=log(value),fill=depth_bin,group=interaction(Time,depth_bin)),position=position_dodge(width=0.8,preserve = "single"),outliers = FALSE)+
    scale_fill_brewer(name = 'depth bins (m)')+
    theme_bw()+
    theme(axis.title.x = element_blank(),
          # axis.text.x = element_text(color = ifelse(levels(data$year) %in% c(2016,2018,2019,2022), "red",
          #                                                                          ifelse(levels(data$year) %in% c(1992,2007:2010,2012), "blue",'black')),size=14),
          plot.title = element_text(size=16),legend.text = element_text(size = 12),  # Increase legend text size
          legend.title = element_text(size = 14), # Increase legend title size
          legend.key.size = unit(1.5, 'lines'),axis.title.y = element_text(size=14),
          # legend.position = 'none'
          )+
    ggtitle(label =sp)+
    scale_x_continuous(limits = c(2001.5,2016.5),breaks=c(2002:2016))+
    ylab(label = 'log(density) (kg/km2)')+
    scale_fill_manual(values = c("0-50" = "#caf0f8",'50-100'='#90e0ef', "100-150" = "#00b4d8",'150-200'='#0077B6','>200'='#023e8a', "warm" = "red", "cold" = "blue"), name = "env . depthbin")

  ########################################################################
  ########################################################################
  
  ########################
  ########################
  # smooth with LOESS
  ########################
  ########################
  
  grid_ebs<-subset(grid.ebs_year,region=='EBSshelf' & Year %in% 1982:2022)
  data_ebs<-cbind(dens_shelf1ebs,'depth'=grid_ebs$DepthGEBCO,'temp'=grid_ebs$Temp)
  
  #merge time region
  data_ebs<-merge(data_ebs,subset(temp_region,region=='EBSshelf'),by.x='Time',by.y='year')
  
  # Define a custom color scale function
  custom_colors <- colorRampPalette(c("blue", "white", "darkred"))
  
  # Plot the original data points
  ggplot(data=subset(data_ebs,value!=0), aes(x = depth, y = value,color=temp.y)) +
    geom_point(alpha=0.5,shape=16) +
    scale_color_gradientn(colors = custom_colors(100),name = 'SBT (°C)')+
    scale_x_continuous(limits=c(0,400))+
    labs(title = paste0(sp, "/ Biomass vs Depth in EBSshelf"), x = "Depth", y = "Biomass")
  
  # Create a list to store the smooth data for each year
  smooth_data_list <- list()
  
  #for selected years
  #for (y in c(2008,2009,2012,2016,2019,2022)) {
  for (y in c(1982:2022)) {
      
    #y<-2008
    
    # Filter data for the current year
    data_ebs_year <- subset(data_ebs,Time == y)
    
    # Apply loess smoothing
    loess_fit <- loess(value ~ depth, data = data_ebs_year, span = 0.5)
    
    # Predict values using the loess model for a smooth curve
    depth_range <- seq(0, max(data_ebs_year$depth), length.out = 100)
    biomass_pred <- predict(loess_fit, newdata = data.frame(depth = depth_range))
    
    # Create a data frame with the predicted values for plotting
    smooth_data <- data.frame(Time = y, depth = depth_range, value = biomass_pred)
    
    # Store the smooth data in the list
    smooth_data_list[[as.character(y)]] <- smooth_data
    
  }
  
  #subset to a few years (representative of cold vs warm)
  #data_ebs1<-subset(data_ebs,Time %in% c(2008,2009,2012,2016,2019,2022))
  
  # Combine all smooth data into one data frame
  smooth_data_all <- do.call(rbind, smooth_data_list)
  
  #merge time region
  smooth_data_all<-merge(smooth_data_all,subset(temp_region,region=='EBSshelf'),by.x='Time',by.y='year')
  
  # Plot the original data points and the smooth curves for each year
  ggplot() +
    geom_point(data=data_ebs, aes(x = depth, y = value, color = temp.y),alpha = 0.2,shape=16) +
    geom_line(data = smooth_data_all, aes(x = depth, y = value, color = temp,group=Time), size = 1) +
    labs(title = paste(sp, "/ Biomass vs Depth with Smooth Curves in EBSshelf for selected years"), x = "Depth", y = "Biomass") +
    theme_minimal() +
    scale_color_gradientn(colors = custom_colors(100),name = 'SBT (°C)')+
    #scale_fill_gradientn(colors = custom_colors(100),name = 'SBT (°C)')+
    scale_x_continuous(limits=c(0,400))+
    theme(text=element_text(size=14),legend.position = c(0.75,0.75))
  
  ########################################################################
  ########################################################################
  
  
}

########################
########################
# metrics using predicted biomass from OM
########################
########################

#df to store results
metrics_df<-data.frame(matrix(NA,nrow = 0,ncol=10))
names(metrics_df)<-c('sp','year','bio','cog_lat','cog_lon','cog_depth','cog_temp','maxdepth','minlat','region')

#loop over common spp in shelf and slope
for (s in sp_shelfslope) {
  
  #s<-sp_shelfslope[1];y<-colnames(dens_slope)[1]
  
  cat(paste('#################### ',s," - ",'\n'))
  
  #EBS+NBS model fit and density
  load(paste0('./shelf EBS NBS VAST/',s,'/fit.RData'))
  fit_shelf<-fit
  bio_shelf<-fit_shelf$Report$Index_gctl[,,,1] #biomass
  dens_shelf<-fit_shelf$Report$D_gct[,1,] #dens
  
  #CHECK units
  #dim(fit_shelf$Report$D_gct[,1,])
  #dim(fit_shelf$Report$Index_gctl[,,,1])
  #x<-fit_shelf$Report$D_gct[,1,] #dens
  #xx<-x*fit_shelf$data_list$a_gl[,1] #bio
  #y<-fit_shelf$Report$Index_gctl[,,,1] #bio
  
  #slope model fit ande density
  load(paste0('./slope EBS VAST/',s,'/fit.RData'))
  fit_slope<-fit
  bio_slope<-fit_slope$Report$Index_gctl[,,,1] #biomass
  bio_slope<-bio_slope*2 #biomass
  dens_slope<-fit_slope$Report$D_gct[,1,] #dens
  dens_slope<-dens_slope*2 #dens
  
  #loop over years
  for (y in colnames(bio_slope)) {
    
    grid.ebs_year1<-subset(grid.ebs_year,Year==y)
    all_bio2<-c(bio_shelf[,y],bio_slope[,y])
    all_dens2<-c(dens_shelf[,y],dens_slope[,y])
    
    #loop over regions
    for (r in c('EBSshelf','NBS','EBS+NBS','EBSslope','all')) {
      
      if (r=='EBSshelf') {
        c<-which(grid.ebs_year1$region=='EBSshelf')
      } else if (r=='NBS') {
        c<-which(grid.ebs_year1$region=='NBS')
      } else if (r == 'EBS+NBS') {
        c <- which(grid.ebs_year1$region %in% c('NBS','EBSshelf'))
      } else if (r=="EBSslope") {
        c<-which(grid.ebs_year1$region=='EBSslope')
      } else if (r=='all') {
        c<-which(!is.na(grid.ebs_year1$region))
      }
      
      #subset grid and dens by grid
      igrid<-grid.ebs_year1[c,]
      all_bio3<-all_bio2[c]
      all_dens3<-all_dens2[c]
      
      #merge dens and grid
      dens_df<-drop_units(cbind(igrid,'dens'=all_dens3))
      
      #only positive
      #dens_df<-dens_df[which(dens_df$dens!=0),]
      
      #percentile
      ipercentile<-t(quantile(dens_df$dens, probs = c(0.10),na.rm=TRUE))[1,]
      
      #data from grid below the percentile
      below_ipercentile<-dens_df[which(dens_df$dens < ipercentile),]
      
      ggplot()+
        geom_density(data=below_ipercentile,aes(x=DepthGEBCO))
      
      
      #max depth and min lat at below density percentile
      maxdepth<-max(below_ipercentile$DepthGEBCO,na.rm = TRUE)
      minlat<-min(below_ipercentile$Lat,na.rm = TRUE)
      
      
      
      #append
      metrics_df<-rbind(metrics_df,
                    data.frame('sp'=s,
                               'year'=y,
                               'bio'=sum(all_bio3,na.rm = TRUE),
                               'cog_lat'=sum(all_dens3*igrid$Lat,na.rm = TRUE)/sum(all_dens3,na.rm = TRUE),
                               'cog_lon'=sum(all_dens3*igrid$Lon,na.rm = TRUE)/sum(all_dens3,na.rm = TRUE),
                               'cog_depth'=sum(all_dens3*igrid$DepthGEBCO,na.rm = TRUE)/sum(all_dens3,na.rm = TRUE),
                               'cog_temp'=sum(all_dens3*igrid$Temp,na.rm = TRUE)/sum(all_dens3,na.rm = TRUE),
                               'maxdepth'=maxdepth,
                               'minlat'=minlat,
                               'region'=r,
                               'percentiles'=data.frame(t(quantile(dens_df$dens[dens_df$dens!=0], probs = c(0.05, 0.10,0.20,0.80, 0.90, 0.95),na.rm=TRUE))))) # Calculate percentiles)))  
      
    }
  }
}

#save table
save(metrics_df,file='./output/metrics_df.RData')
#load(file='./output/metrics_df.RData')

########################
########################
#  depth and lat cog of temp
########################
########################

#df to store results
metrics_base<-data.frame(matrix(NA,nrow = 0,ncol=3))
names(metrics_base)<-c('year','cog_lat','cog_depth')

for (y in 1982:2022) {
  
  #y<-1982
  
  igrid<-subset(grid.ebs_year,Year==y)
  
  metrics_base<-rbind(metrics_base,
                    data.frame('year'=y,
                               'cog_lat'=sum(igrid$Temp*igrid$Lat,na.rm = TRUE)/sum(igrid$Temp,na.rm = TRUE),
                               'cog_depth'=sum(igrid$Temp*igrid$depth_m,na.rm = TRUE)/sum(igrid$Temp,na.rm = TRUE))) # Calculate percentiles)))  
  
  
}

#save table
save(metrics_base,file='./output/metrics_base.RData')
load(file='./output/metrics_base.RData')

#temp region
metrics_base1<-merge(metrics_base,subset(temp_region,region=='all'),by=c('year'))

p<-
  ggplot()+
  #geom_point(data=metrics_df2,aes(x=cog_depth,y=cog_lat,fill=temp),shape=21)+
  geom_text(data=subset(metrics_base1,year %in% bold_years),aes(x=cog_depth,y=cog_lat),label='*',color='black',fontface='bold',nudge_x = 0.95,nudge_y = 0.03)+
  geom_text(data=subset(metrics_base1,year %in% 2002:2016),aes(x=cog_depth,y=cog_lat,color=temp,label=year),fontface='bold')+
  labs(x='bathymetrical COG relative to SBT',y='latitudinal COG relative to SBT')+
  theme_bw()+
  theme(aspect.ratio = 1)+
  scale_color_gradientn(colors = custom_colors(100),name = 'SBT (°C)')

agg_png(paste0('./figures slope/base_cog.png'), width = 4, height = 4, units = "in", res = 300)
p
dev.off()

########################
########################
# cog lat vs cog depth
########################
########################

#temp region
metrics_df1<-merge(metrics_df,temp_region,by=c('year','region'))

#drop units
metrics_df1<-drop_units(metrics_df1)

# Define a custom color scale function
custom_colors <- colorRampPalette(c("blue", "white", "darkred"))

#loop over common spp in shelf and slope
for (s in sp_shelfslope) {

  s<-sp_shelfslope[1]
  
  #subset
  metrics_df2<-subset(metrics_df1,sp==s & region=='all')
  
  p1a<-
  ggplot()+
    #geom_point(data=metrics_df2,aes(x=cog_depth,y=cog_lat,fill=temp),shape=21)+
    geom_text(data=subset(metrics_df2,year %in% bold_years),aes(x=cog_depth,y=cog_lat),label='*',color='black',fontface='bold',nudge_x = 0.95,nudge_y = 0.02)+
    geom_text(data=metrics_df2,aes(x=cog_depth,y=cog_lat,color=temp,label=year),fontface='bold')+
    labs(x='bathymetrical COG',y='latitudinal COG')+
    theme_bw()+
    theme(legend.position = 'none',aspect.ratio = 1)+
    scale_color_gradientn(colors = custom_colors(100),name = 'SBT (°C)')
  
  p1b<-
  ggplot()+
    #geom_point(data=metrics_df2,aes(x=cog_depth,y=cog_lat,fill=temp),shape=21)+
    geom_text(data=subset(metrics_df2,year %in% bold_years),aes(x=cog_depth,y=cog_lat),label='*',color='black',fontface='bold',nudge_x = 0.95,nudge_y = 0.02)+
    geom_text(data=metrics_df2,aes(x=cog_depth,y=cog_lat,color=log(bio),label=year),fontface='bold')+
    labs(x='bathymetrical COG',y='latitudinal COG')+
    theme_bw()+
    theme(legend.position = 'none',aspect.ratio = 1)+
    scale_color_viridis_c(name = 'log(bio)')
  
  p2a<-
  ggplot()+
    #geom_point(data=metrics_df2,aes(x=cog_depth,y=cog_lat,fill=temp),shape=21)+
    geom_text(data=subset(metrics_df2,year %in% bold_years),aes(x=maxdepth,y=minlat),label='*',color='black',fontface='bold',nudge_x = 0.95,nudge_y = 0.02)+
    geom_text(data=metrics_df2,aes(x=maxdepth,y=minlat,color=temp,label=year),fontface='bold')+
    labs(x='max depth below P10 bio',y='min lat  below P10 bio')+
    theme_bw()+
    theme(aspect.ratio = 1)+
    scale_color_gradientn(colors = custom_colors(100),name = 'SBT (°C)',guide = guide_colorbar(frame.colour = "black", ticks.colour = "black"))   
  
  p2b<-
  ggplot()+
    #geom_point(data=metrics_df2,aes(x=cog_depth,y=cog_lat,fill=temp),shape=21)+
    geom_text(data=subset(metrics_df2,year %in% bold_years),aes(x=maxdepth,y=minlat),label='*',color='black',fontface='bold',nudge_x = 0.95,nudge_y = 0.02)+
    geom_text(data=metrics_df2,aes(x=maxdepth,y=minlat,color=log(bio),label=year),fontface='bold')+
    labs(x='max depth  below P10 bio',y='min lat below P10 bio')+
    theme_bw()+
    theme(aspect.ratio = 1)+
    scale_color_viridis_c(name = 'log(bio)',guide = guide_colorbar(frame.colour = "black", ticks.colour = "black"))
  
  p3a<-
    ggplot()+
    #geom_point(data=metrics_df2,aes(x=cog_depth,y=cog_lat,fill=temp),shape=21)+
    geom_text(data=subset(metrics_df2,year %in% bold_years),aes(x=cog_depth,y=cog_temp),label='*',color='black',fontface='bold',nudge_x = 0.95,nudge_y = 0.02)+
    geom_text(data=metrics_df2,aes(x=cog_depth,y=cog_temp,color=temp,label=year),fontface='bold')+
    labs(x='bathymetrical COG',y='temp COG')+
    theme_bw()+
    theme(legend.position = 'none',aspect.ratio = 1)+
    scale_color_gradientn(colors = custom_colors(100),name = 'SBT (°C)')
  
  p3b<-
    ggplot()+
    #geom_point(data=metrics_df2,aes(x=cog_depth,y=cog_lat,fill=temp),shape=21)+
    geom_text(data=subset(metrics_df2,year %in% bold_years),aes(x=cog_depth,y=cog_temp),label='*',color='black',fontface='bold',nudge_x = 0.95,nudge_y = 0.02)+
    geom_text(data=metrics_df2,aes(x=cog_depth,y=cog_temp,color=log(bio),label=year),fontface='bold')+
    labs(x='bathymetrical COG',y='temp COG')+
    theme_bw()+
    theme(legend.position = 'none',aspect.ratio = 1)+
    scale_color_viridis_c(name = 'log(bio)')

  
  #title
  title <- ggdraw() + draw_label(s, fontface='bold')
  
  #plots
  plots<-
    plot_grid(p1a,p3a,p2a,p1b,p3b,p2b,nrow = 2,rel_widths = c(0.9,0.9,1.19))
  
  #save  plot
  agg_png(paste0('./figures slope/shift_',s,'_c2','.png'), width = 9, height = 6, units = "in", res = 300)
  print(
  plot_grid(title,plots,ncol=1, rel_heights=c(0.05, 1))
  )
  dev.off()
  
}

########################
########################
# range edges gif
########################
########################

# Initialize variables
output_dir <- "frames/"
dir.create(output_dir, showWarnings = FALSE)  # Create directory for frames

# Loop over species
for (s in sp_shelfslope) {
  
  s<-sp_shelfslope[1]
  
  # Load shelf model
  load(paste0('./shelf EBS NBS VAST/', s, '/fit.RData'))
  fit_shelf <- fit
  bio_shelf <- fit_shelf$Report$Index_gctl[,,,1]
  
  # Load slope model
  load(paste0('./slope EBS VAST/', s, '/fit_st.RData'))
  fit_slope <- fit
  bio_slope <- fit_slope$Report$Index_gctl[,,,1]
  
  #to save area output
  area_df<-data.frame(matrix(NA,nrow = 0,ncol=5))
  colnames(area_df)<-c('region','abundance','mean_dens','area_occupied','year')
  
  # Loop over years
  for (y in as.factor(2002:2016)) {
    
    cat(paste('#################### ', s, " - ", y, '\n'))
    
    grid.ebs_year1 <- subset(grid.ebs_year, Year == y)
    
    all_bio2 <- c(bio_shelf[, y], bio_slope[, y])
    all_bio3 <- cbind(all_bio2, grid.ebs_year1)
    all_bio3$dens <- all_bio3$all_bio2 / all_bio3$Area_in_survey_km2
    all_bio4 <- subset(all_bio3, region %in% c('EBSshelf', 'EBSslope', 'NBS') & DepthGEBCO > 0)
    
    ########################
    ########################
    # effective area occuppied
    ########################
    ########################
    
    #total abundance bt is equal to average density mt 
    #ht=bt+mt
    
    #mean density
    all_meandens<-mean(all_bio3$dens)
    region_meandens<-aggregate(dens ~ region,all_bio3,FUN=mean,na.rm=TRUE)
    region_meandens<-rbind(region_meandens,c('all',all_meandens))
    #total abundance
    all_abu<-sum(all_bio3$all_bio2)
    region_abu<-aggregate(all_bio2 ~ region,all_bio3,FUN=sum,na.rm=TRUE)
    region_abu<-rbind(region_abu,c('all',all_abu))
    #append df
    area_df1<-
    data.frame('region'=region_meandens$region,
               'abundance'=region_abu$all_bio2,
               'mean_dens'=region_meandens$dens,
               'area_occupied'=as.numeric(region_abu$all_bio2)/as.numeric(region_meandens$dens),
               'year'=y)
    
    area_df<-rbind(area_df,area_df1)
    
    ########################
    ########################
    # range edges
    ########################
    ########################
    
    p5 <- quantile(all_bio4$dens, c(0.10))
    p95 <- quantile(all_bio4$dens, c(0.90))
    
    all_bio4$quantile <- ifelse(all_bio4$dens < p95 & all_bio4$dens > p5, 1, 0)
    
    coordinates(all_bio4) <- ~ Lon + Lat
    crs(all_bio4) <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
    
    all_bio5 <- spTransform(all_bio4, '+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs')
    all_bio6 <- data.frame(all_bio5)
    
    xycells <- as.integer(sqrt(dim(all_bio5)[1]))
    
    r1 <- raster(ext = extent(all_bio5), ncol = xycells, nrow = xycells)
    
    r2 <- rasterize(all_bio5, r1, field = c('quantile', 'Temp', 'all'))
    crs(r2) <- CRS('+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs')
    r2$quantile[r2$quantile == 0] <- NA
    r2$Temp[r2$Temp > 2] <- NA
    r2[r2$Temp < 2] <- 1
    
    r3 <- as.data.frame(r2, xy = TRUE)
    r3 <- r3[complete.cases(r3$quantile),]
    r3$quantile <- as.factor(r3$quantile)
    
    r4a <- rasterToPolygons(r2$quantile, dissolve = TRUE, digits = 1)
    r4b <- rasterToPolygons(r2$Temp, dissolve = TRUE, digits = 1)
    #r4c <- rasterToPolygons(r2$all, dissolve = TRUE, digits = 1)
    
    get_temp_color <- function(itemp, temp_min, temp_max) {
      gradient <- colorRampPalette(c("blue", "white", "darkred"))(100)
      temp_scaled <- scales::rescale(itemp, to = c(1, 100), from = c(temp_min, temp_max))
      return(gradient[round(temp_scaled)])
    }
    
    alltemp_region <- subset(temp_region, region == 'all' & year %in% 2002:2016)
    temp_min <- min(alltemp_region$temp)
    temp_max <- max(alltemp_region$temp)
    
    itemp_df <- subset(alltemp_region, year == y)
    itemp <- itemp_df$temp
    
    # Generate plot
    p <- 
      ggplot() +
      geom_polygon(data = r4b, aes(x = long, y = lat, group = group, color = 'coldpool'), fill = 'blue', alpha = 0.7) +
      geom_tile(data = r3, aes(x = x, y = y, fill = quantile), alpha = 0.4) +
      scale_fill_manual(values = 'grey30', labels = 'species range') +
      geom_polygon(data = r4a, aes(x = long, y = lat, group = group), fill = 'transparent', colour = "grey50") +
      scale_color_manual(values = c('cold pool' = 'transparent')) +
      geom_polygon(data = ak_sppoly, aes(x = long, y = lat, group = group), color = 'transparent', linewidth = 0.2, fill = 'white') +
      geom_polygon(data = NBS_sh, aes(x = long, y = lat, group = group), fill = NA, col = 'black') +
      geom_polygon(data = EBSshelf_sh, aes(x = long, y = lat, group = group), fill = NA, col = 'black') +
      geom_polygon(data = EBSslope_sh, aes(x = long, y = lat, group = group), fill = NA, col = 'black') +
      theme_bw() +
      labs(title = s) +
      geom_shadowtext(aes(x = -1200000, y = 1750000, label = y),color = get_temp_color(itemp, temp_min, temp_max), size = 14, fontface = "bold", bg.colour = "black", bg.r = .05) +
      #geom_label(aes(x = -1200000, y = 1750000, label = y), size = 10,color='black', fill = get_temp_color(itemp, temp_min, temp_max),) +
      #geom_label(aes(x = -1200000, y = 1750000, label = y), size = 10,fill='black', color = get_temp_color(itemp, temp_min, temp_max),) +
      #annotate(geom="text",x = -1000000, y = 1900000, label = s, size = 10) +
      coord_sf(crs = '+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs',
               xlim = c(panel_extent$x[1] + 200000, panel_extent$x[2] + 30000),
               ylim = c(panel_extent$y[1] - 100000, panel_extent$y[2] - 200000),
               label_axes = "-NE-") +
      theme(axis.title = element_blank(),
            legend.title = element_blank(),
            plot.margin = unit(c(0, 0, 0, 0), units = 'in'),
            axis.text = element_blank(),
            legend.position = 'none',
            panel.grid = element_blank(),
            plot.title = element_text(hjust = 0.5,vjust = -12,size = 18),
            panel.border = element_blank(),
            axis.ticks = element_blank())
    
    # Save each plot as a frame
    frame_path <- paste0(output_dir, s, "_", y, ".png")
    ggsave(frame_path, plot = p, width = 8, height = 6)
    
    
  }
  
  # Combine frames into a GIF
  frame_files <- list.files(output_dir, pattern = s, full.names = TRUE)
  gif_path <- paste0("./figures slope/range_edge_",s,".gif")
  images <- lapply(frame_files, image_read)
  image_join(images) %>%
    image_animate(fps = 1) %>%
    image_write(gif_path)
  
  # Define a data frame for the rectangles
  rect_data <- data.frame(
    xmin = c(2002, 2005.5, 2013.5),
    xmax = c(2005.5, 2013.5, 2016),
    fill = factor(c("warm", "cold", "warm"))
  )
  
  #area df changes
  ggplot()+
    #geom_point(data=metrics_df2,aes(x=cog_depth,y=cog_lat,fill=temp),shape=21)+
    geom_rect(data = rect_data, aes(ymin = -Inf, ymax = Inf, xmin = xmin, xmax = xmax, fill = fill), alpha = 0.3) +
    geom_line(data=area_df,aes(x=as.numeric(year),y=log(area_occupied),group=region),color='black')+
    labs(x='',y='effective area occupied',title=s)+
    theme_bw()+
    facet_wrap(~region,scales = 'free_y',nrow = 1)+
    theme(legend.position = 'none',aspect.ratio = 1)
    #scale_color_viridis_c(name = 'log(bio)')
  
  
}



