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
pack_cran<-c('cowplot','ggspatial','raster','rasterVis','rgeos','scales','rnaturalearth','grid','ggplot2','lubridate','ragg','rgdal','dplyr','sp')

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

#get detph in grid
grid.ebs_year1<-unique(grid.ebs_year[,c("DepthGEBCO","region","Area_in_survey_km2","Lat","Lon")])
#grid.ebs_year1s<-unique(grid.ebs_year[,c('Stratum',"DepthGEBCO","region","Area_in_survey_km2","Lat","Lon")])

#add depths regions
region_depth<-ifelse(grid.ebs_year1$DepthGEBCO<100,'0_100',
                                    ifelse(grid.ebs_year1$DepthGEBCO>=100 & grid.ebs_year1$DepthGEBCO<=200,'100_200',
                                           ifelse(grid.ebs_year1$DepthGEBCO>=200 & grid.ebs_year1$DepthGEBCO<=300,'200_300',
                                                  ifelse(grid.ebs_year1$DepthGEBCO>=300 & grid.ebs_year1$DepthGEBCO<=400,'300_400',NA))))

region_depth<-ifelse(grid.ebs_year1$region=='NBS',NA,region_depth)

#cbind
grid.ebs_year1<-cbind(grid.ebs_year1,region_depth)

#plot
ggplot()+
  geom_point(data=grid.ebs_year1,aes(x=Lon,y=Lat,color=region_depth))

#average depth by region and year
depth_region<-
  aggregate(DepthGEBCO ~ region + region,grid.ebs_year,FUN=mean)
grid.ebs_year_d1<-subset(grid.ebs_year,region %in% c('EBSshelf','NBS'))
depth_region<-
rbind(depth_region,
      c('EBS+NBS',mean(grid.ebs_year_d1$DepthGEBCO)),
      c('all',mean(grid.ebs_year$DepthGEBCO)))


#average SBT by region and year
temp_region<-
  aggregate(Temp ~ Year + region,grid.ebs_year,FUN=mean)
names(temp_region)<-c('year','region','temp')

ebsnbstemp<-
  aggregate(Temp ~ Year ,subset(grid.ebs_year,region %in% c('EBSshelf','NBS')),FUN=mean)
ebsnbstemp<-data.frame('year'=ebsnbstemp$Year,
                       'region'='EBS+NBS',
                       'temp'=ebsnbstemp$Temp)

#get average SBT for all region (EBSshelf+NBS+EBSslope)
alltemp<-aggregate(temp ~ year, temp_region,FUN=mean)
alltemp$region<-'all'
alltemp<-alltemp[,c("year","region","temp")]

#rbind
temp_region<-rbind(temp_region,ebsnbstemp,alltemp)

#plot
p<-
ggplot()+
  geom_rect(aes(xmin=2002,xmax=2016,ymin=-Inf,ymax=Inf),alpha=0.1)+
  geom_point(data=temp_region,aes(x=year,y=temp,color=region),alpha=0.7)+
  scale_color_manual(values=c('EBSslope'='#B4AF46','EBSshelf'='#B4464B','NBS'='#4682B4','EBS+NBS'='#b649b1','all'='black'),
                     breaks = c('all','EBS+NBS','EBSshelf','NBS','EBSslope'),
                     labels=c('ALL','EBS+NBS','EBSshelf','NBS','EBSslope'))+
  theme_bw()+
  labs(x='',y='average SBT')+
  geom_vline(xintercept = c(2003:2005,2015,2016),linetype='dashed',color='red',alpha=0.5)+
  geom_vline(xintercept = c(2008,2009,2012),linetype='dashed',color='blue',alpha=0.5)

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
temp_region11<-subset(temp_region1,region %in% c('EBSshelf','EBS+NBS','NBS'))

#cold/warm years
cyrs<-unique(temp_region11[which(temp_region11$typeyear %in% c('cold')),'year'])
wyrs<-unique(temp_region11[which(temp_region11$typeyear %in% c('warm')),'year'])

#plot
p<-
ggplot()+
  geom_rect(aes(xmin=1982,xmax=2022,ymin=-1,ymax=+1),linetype='dashed',alpha=0.2)+
  geom_vline(xintercept = cyrs,linetype='dashed',color='blue',alpha=0.5)+
  geom_vline(xintercept = wyrs,linetype='dashed',color='red',alpha=0.5)+
  geom_hline(yintercept = 0,alpha=0.5,linetype='dashed')+
  geom_point(data=temp_region1,aes(x=year,y=scaledtemp,color=region),alpha=0.7)+
  scale_color_manual(values=c('EBSslope'='#B4AF46','EBSshelf'='#B4464B','NBS'='#4682B4','EBS+NBS'='#b649b1','all'='black'),
                     breaks = c('all','EBS+NBS','EBSshelf','NBS','EBSslope'),
                     labels=c('ALL','EBS+NBS','EBSshelf','NBS','EBSslope'))+
  theme_bw()+
  scale_x_continuous(limits = c(1982,2022))+
  labs(x='',y='average scaled SBT')

  #geom_vline(xintercept = c(2008,2009,2012),linetype='dashed',color='blue',alpha=0.5)

#save  plot
agg_png(paste0('./figures slope/scaledSBT_regions.png'), width = 7, height = 4, units = "in", res = 300)
print(p)
dev.off()

temp_region1[which(temp_region1$typeyear %in% c('warm','cold')),'year']

########################
########################
# using predicted biomass from OM
########################
########################

#df to store results
cog_df<-data.frame(matrix(NA,nrow = 0,ncol=7))
names(cog_df)<-c('sp','year','bio','cog_lat','cog_lon','cog_depth','region')

#loop over common spp in shelf and slope
for (s in sp_shelfslope) {
  
  #s<-sp_shelfslope[1];y<-colnames(dens_slope)[1]
  
  cat(paste('#################### ',s," - ",'\n'))
  
  #EBS+NBS model fit and density
  load(paste0('./shelf EBS NBS VAST/',s,'/fit.RData'))
  fit_shelf<-fit
  bio_shelf<-fit_shelf$Report$Index_gctl[,,,1] #biomass
  
  #CHECK units
  dim(fit_shelf$Report$D_gct[,1,])
  dim(fit_shelf$Report$Index_gctl[,,,1])
  x<-fit_shelf$Report$D_gct[,1,] #dens
  xx<-x*fit_shelf$data_list$a_gl[,1] #bio
  y<-fit_shelf$Report$Index_gctl[,,,1] #bio
  
  #slope model fit ande density
  load(paste0('./slope EBS VAST/',s,'/fit.RData'))
  fit_slope<-fit
  bio_slope<-fit_slope$Report$Index_gctl[,,,1] #biomass

  #loop over years
  for (y in colnames(bio_slope)) {
    
    all_bio2<-c(bio_shelf[,y],bio_slope[,y])
    
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
    
      #append
      cog_df<-rbind(cog_df,
                    data.frame('sp'=s,
                               'year'=y,
                               'bio'=sum(all_bio3),
                               'cog_lat'=sum(all_bio3*igrid$Lat)/sum(all_bio3),
                               'cog_lon'=sum(all_bio3*igrid$Lon)/sum(all_bio3),
                               'cog_depth'=sum(all_bio3*igrid$DepthGEBCO)/sum(all_bio3),
                               'region'=r))  
          
    }
  }
}

#save table
save(cog_df,file='./output/cog_OM.RData')
#load(file='./output/cog_OM.RData')

###################
# merge 
###################

cog_df1<-merge(cog_df,temp_region,by=c('year','region'))

#warm<-c(2002:2005,2014:2016)
#cold<-c(2006:2013)
cog_df1$coloryear<-ifelse(cog_df1$year %in% c(2002:2005,2014:2016),'red','blue')

###################
# plot region temp vs depth cog over species and regions 
###################

# Define a function to scale values
scale_value <- function(x) {
  (x - mean(x)) / sd(x)
}


# Add a scaled column using group_by and mutate
cog_df11 <- cog_df1 %>%
  group_by(sp, region) %>%
  mutate(scale_bio = scale_value(bio))

p<-
ggplot() +
  geom_point(data = cog_df11, aes(x = cog_depth, y = temp, color = as.factor(year),size=scale_bio),alpha=0.6) +
  facet_wrap(region ~ sp, scales = 'free', nrow = 5) +
  labs(x = 'bathymetric center of gravity', y = 'mean region SBT') +
  scale_color_manual(values = unique(cog_df1[,c("coloryear","year")])[,'coloryear'], 
                     breaks = unique(cog_df1[,c("coloryear","year")])[,'year'],  # Specify the breaks based on unique year values
                     labels = unique(cog_df1[,c("coloryear","year")])[,'year'],  # Use unique year values as labels
                     name = 'year') +
  #geom_vline(data=depth_region,aes(xintercept=as.numeric(DepthGEBCO)),linetype='dotted')+
  theme_bw()+
  theme(legend.position = 'none')+
  #geom_smooth(data = cog_df1, aes(x = depth, y = temp), method = "lm", se = FALSE, color = 'black',alpha=0.9)
  geom_line(data = cog_df1, aes(x = cog_depth, y = temp),stat="smooth",method = "lm", 
            linetype ="dashed",
            color='black')+
  scale_size() #limits = c(100000,1000000000000),

#save  plot
agg_png(paste0('./figures slope/regionSBT_cogdepth.png'), width = 18, height =8, units = "in", res = 300)
print(p)
dev.off()

########################
########################
# using predicted biomass from OM + region
########################
########################

#add depths regions
region_depth1<-ifelse(grid.ebs_year1$DepthGEBCO<50,'coastal',
                      ifelse(grid.ebs_year1$DepthGEBCO>=50 & grid.ebs_year1$DepthGEBCO<=100,'middle',
                             ifelse(grid.ebs_year1$DepthGEBCO>=100 & grid.ebs_year1$DepthGEBCO<=200,'outer',NA)))

region_depth1<-ifelse(grid.ebs_year1$region=='NBS',NA,region_depth1)

#cbind
grid.ebs_year1<-cbind(grid.ebs_year1,region_depth1)

#df to store results
cog_df<-data.frame(matrix(NA,nrow = 0,ncol=7))
names(cog_df)<-c('sp','year','bio','cog_lat','cog_lon','cog_depth','region')

#loop over common spp in shelf and slope
for (s in sp_shelfslope) {
  
  #s<-sp_shelfslope[1];y<-colnames(dens_slope)[1]
  
  cat(paste('#################### ',s," - ",'\n'))
  
  #EBS+NBS model fit and density
  load(paste0('./shelf EBS NBS VAST/',s,'/fit.RData'))
  fit_shelf<-fit
  bio_shelf<-fit_shelf$Report$Index_gctl[,,,1] #biomass
  
  #CHECK units
  dim(fit_shelf$Report$D_gct[,1,])
  dim(fit_shelf$Report$Index_gctl[,,,1])
  x<-fit_shelf$Report$D_gct[,1,] #dens
  xx<-x*fit_shelf$data_list$a_gl[,1] #bio
  y<-fit_shelf$Report$Index_gctl[,,,1] #bio
  
  # #slope model fit ande density
  # load(paste0('./slope EBS VAST/',s,'/fit.RData'))
  # fit_slope<-fit
  # bio_slope<-fit_slope$Report$Index_gctl[,,,1] #biomass
  # 
  #loop over years
  for (y in colnames(bio_shelf)) {
    
    #all_bio2<-c(bio_shelf[,y],bio_slope[,y])
    all_bio2<-c(bio_shelf[,y])
    
    #loop over regions
    for (r in c('coastal','middle','outer')) {
      
      if (r=='coastal') {
        c<-which(grid.ebs_year1$region_depth1=='coastal')
      } else if (r=='middle') {
        c<-which(grid.ebs_year1$region_depth1=='middle')
      } else if (r == 'outer') {
        c <- which(grid.ebs_year1$region_depth1=='outer')
      } 
      
      #subset grid and dens by grid
      igrid<-grid.ebs_year1[c,]
      all_bio3<-all_bio2[c]
      
      #append
      cog_df<-rbind(cog_df,
                    data.frame('sp'=s,
                               'year'=y,
                               'bio'=sum(all_bio3,na.rm = TRUE),
                               'cog_lat'=sum(all_bio3*igrid$Lat,na.rm = TRUE)/sum(all_bio3,na.rm = TRUE),
                               'cog_lon'=sum(all_bio3*igrid$Lon,na.rm = TRUE)/sum(all_bio3,na.rm = TRUE),
                               'cog_depth'=sum(all_bio3*igrid$DepthGEBCO,na.rm = TRUE)/sum(all_bio3,na.rm = TRUE),
                               'region'=r))  
      
    }
  }
}

#save table
save(cog_df,file='./output/cog_OM_ebs.RData')
#load(file='./output/cog_OM_ebs.RData')

###################
# merge 
###################

ebs_grid.ebs_year<-subset(grid.ebs_year,region=='EBSshelf')
ebs_grid.ebs_year$region_depth<-ifelse(ebs_grid.ebs_year$Stratum %in% c('50','61','62','90'),'outer',
                                        ifelse(ebs_grid.ebs_year$Stratum %in% c('10','20'),'coastal','middle'))

ggplot()+
  geom_point(data=subset(ebs_grid.ebs_year,Year=='1982'),aes(x=Lon,y=Lat,color=region_depth))
  # scale_color_manual(values=c('50'='blue','61'='blue','62'='blue','90'='blue',
  #                             '31'='green','32'='green','41'='green','42'='green','43'='green','82'='green',
  #                             '10'='red','20'='red'))

#average SBT by region and year
ebs_temp_region<-
  aggregate(Temp ~ Year + region_depth,ebs_grid.ebs_year,FUN=mean)
names(ebs_temp_region)<-c('year','region','temp')

#average SBT by region and year
ebs_temp_region1<-
  aggregate(Temp ~ Year ,ebs_grid.ebs_year,FUN=mean)
names(ebs_temp_region1)<-c('year','temp')

#cog_df1<-merge(cog_df,temp_region,by=c('year','region'))
cog_df1<-merge(cog_df,ebs_temp_region1,by=c('year'))

#warm<-c(2002:2005,2014:2016)
#cold<-c(2006:2013)
#cog_df1$coloryear<-ifelse(cog_df1$year %in% c(2002:2005,2014:2016),'red','blue')
cog_df1$coloryear<-ifelse(cog_df1$year %in% c(2019,2018,2016),'red',
                          ifelse(cog_df1$year %in% c(1992,1999,2007:2010,2011),'blue','grey'))

###################
# plot region temp vs depth cog over species and regions 
###################

# Define a function to scale values
scale_value <- function(x) {
  (x - mean(x)) / sd(x)
}


# Add a scaled column using group_by and mutate
cog_df11 <- cog_df1 %>%
  group_by(sp, region) %>%
  mutate(scale_bio = scale_value(bio))

p<-
  ggplot() +
  geom_point(data = cog_df11, aes(x = cog_depth, y = temp, color = as.factor(year),size=scale_bio),alpha=0.6) +
  facet_wrap(region ~ sp, scales = 'free', nrow = 3) +
  labs(x = 'bathymetric center of gravity', y = 'mean region SBT') +
  scale_color_manual(values = unique(cog_df1[,c("coloryear","year")])[,'coloryear'], 
                     breaks = unique(cog_df1[,c("coloryear","year")])[,'year'],  # Specify the breaks based on unique year values
                     labels = unique(cog_df1[,c("coloryear","year")])[,'year'],  # Use unique year values as labels
                     name = 'year') +
  #geom_vline(data=depth_region,aes(xintercept=as.numeric(DepthGEBCO)),linetype='dotted')+
  theme_bw()+
  theme(legend.position = 'none')+
  #geom_smooth(data = cog_df1, aes(x = depth, y = temp), method = "lm", se = FALSE, color = 'black',alpha=0.9)
  geom_line(data = cog_df1, aes(x = cog_depth, y = temp),stat="smooth",method = "lm", 
            linetype ="dashed",
            color='black')+
  scale_size() #limits = c(100000,1000000000000),

#save  plot
agg_png(paste0('./figures slope/regionSBT_cogdepthregion.png'), width = 18, height =8, units = "in", res = 300)
print(p)
dev.off()



########################
########################
# using simulated densities
########################
########################

#load simulated densities
load('/Users/daniel/Work/Adapting Monitoring to a Changing Seascape/output/species/ms_sim_dens.RData') #index_hist
sim_dens1_shelf<-sim_dens1
load('/Users/daniel/Work/Adapting Monitoring to a Changing Seascape/output/species/ms_sim_dens_slope.RData') #index_hist
sim_dens1_slope<-sim_dens1
dim(sim_dens1_shelf);dim(sim_dens1_slope)

#to store results
dens_df<-data.frame(matrix(NA,nrow = 0,ncol=9))
names(dens_df)<-c('sp','sim','year','dens','depth','cog_lat','cog_lon','cog_depth','region')

# #loop over common spp in shelf and slope
# for (s in sp_shelfslope) {
#   
#   #s<-sp_shelfslope[1];si<-dimnames(sim_dens1_shelf)[[4]][1];y<-dimnames(sim_dens1_slope)[[3]][1];r<-'EBSshelf'
#   
#   #loop over simulations
#   for (si in dimnames(sim_dens1_shelf)[[4]]) {
#     
#     cat(paste('#################### ',s," - ",si,'\n'))
#     
#     #loop over years
#     for (y in dimnames(sim_dens1_slope)[[3]]) {
#       
#       sim_dens2<-c(sim_dens1_shelf[,s,y,si],sim_dens1_slope[,s,y,si])
#       
#       #regions
#       for (r in c('EBSshelf','NBS','EBSslope','all')) {
#         
#         #rows by region
#         if (r=='EBSshelf') {
#           c<-which(grid.ebs_year1$region=='EBSshelf')
#         } else if (r=='NBS') {
#           c<-which(grid.ebs_year1$region=='NBS')
#         } else if (r=="EBSslope") {
#           c<-which(grid.ebs_year1$region=='EBSslope')
#         } else if (r=='all') {
#           c<-which(!is.na(grid.ebs_year1$region))
#         }
#         
#         #grid
#         igrid<-grid.ebs_year1[c,]
#         sim_dens3<-sim_dens2[c]
#         
#         #store densities
#         dens_df<-rbind(dens_df,
#                        data.frame('sp'=s,
#                                   'sim'=si,
#                                   'year'=y,
#                                   'dens'=sim_dens3,
#                                   'depth'=grid.ebs_year$DepthGEBCO[c],                                
#                                   'cog_lat'=sum(sim_dens3*igrid$Lat)/sum(sim_dens3),
#                                   'cog_lon'=sum(sim_dens3*igrid$Lon)/sum(sim_dens3),
#                                   'cog_depth'=sum(sim_dens3*igrid$DepthGEBCO)/sum(sim_dens3),
#                                   'region'=grid.ebs_year1$region[c])) 
#       }
#     }
#   }
# }

# Indexing vectors outside loops
sims <- dimnames(sim_dens1_shelf)[[4]]
yrs <- dimnames(sim_dens1_slope)[[3]]
yrs<-as.character(1982:2022)

# Pre-allocate memory for cog_df
num_rows <- length(sp_shelfslope)*length(sims)*
  ((length(dimnames(sim_dens1_slope)[[3]])*3)+(length(setdiff(yrs,dimnames(sim_dens1_slope)[[3]]))*5))

cog_df <- data.frame('sp' = character(num_rows),
                      'sim' = character(num_rows),
                      'year' = character(num_rows),
                      'cog_lat' = numeric(num_rows),
                      'cog_lon' = numeric(num_rows),
                      'cog_depth' = numeric(num_rows),
                      'region' = character(num_rows))

#sta<-Sys.time()

# Initialize start indexInitialize start index
start <- 1

for (s in sp_shelfslope) {
  
  for (si in sims) {
    cat(paste('#################### ', s, " - ", si, '\n'))
    
    for (y in yrs) {
      
      if (y %in% dimnames(sim_dens1_slope)[[3]]) {
        
        sim_dens2 <- c(sim_dens1_shelf[, s, y, si], sim_dens1_slope[, s, y, si])
        #sim_bio2<- sim_dens2*grid.ebs_year1$Area_in_survey_km2
        regions<-c('EBSshelf', 'NBS','EBS+NBS', 'EBSslope', 'all')
        
      } else {
        
        sim_dens2 <- c(sim_dens1_shelf[, s, y, si])
        #sim_bio2<- sim_dens2*subset(grid.ebs_year1,region %in% c('EBSshelf','NBS'))[,'Area_in_survey_km2']
        regions<-c('EBSshelf', 'NBS', 'EBS+NBS')
        
      }
        
        
        for (r in regions) {
          
          if (r == 'EBSshelf') {
            c <- which(grid.ebs_year1$region == 'EBSshelf')
          } else if (r == 'NBS') {
            c <- which(grid.ebs_year1$region == 'NBS')
          } else if (r == 'EBS+NBS') {
            c <- which(grid.ebs_year1$region %in% c('NBS','EBSshelf'))
          } else if (r == "EBSslope") {
            c <- which(grid.ebs_year1$region == 'EBSslope')
          } else if (r == 'all') {
            c <- which(!is.na(grid.ebs_year1$region))
          }
          
          igrid <- grid.ebs_year1[c,]
          #sim_bio3 <- sim_bio2[c]
          sim_dens3 <- sim_dens2[c]
          
          cog_df[start:(start+1),] <- data.frame('sp'=s,
                                                            'sim'=si,
                                                            'year'=y,
                                                            'cog_lat'=sum(sim_dens3*igrid$Lat)/sum(sim_dens3),
                                                            'cog_lon'=sum(sim_dens3*igrid$Lon)/sum(sim_dens3),
                                                            'cog_depth'=sum(sim_dens3*igrid$DepthGEBCO)/sum(sim_dens3),
                                                            'region'=r)
        
        
        start <- start + 1
      }
    }
  }
}

#end<-Sys.time()

#save results
save(cog_df,file='./output/cog_simdens.RData')
load(file='./output/cog_simdens.RData')

###################
# merge 
###################

cog_df1<-merge(cog_df,temp_region1,by=c('year','region'))

#warm<-c(2002:2005,2014:2016)
#cold<-c(2006:2013)
#cog_df1$coloryear<-ifelse(cog_df1$year %in% c(2002:2005,2014:2016),'red','blue')

###################
# plot region temp vs depth cog over species and regions 
###################

cog_df1$region<-factor(cog_df1$region,levels = c('all','EBS+NBS','EBSshelf','NBS','EBSslope'))

library(ggpubr)

rsquared_data <- cog_df1 %>%
  na.omit() %>%
  group_by(region, sp) %>%
  do(data.frame(rsquared = summary(lm(temp ~ cog_depth + cog_depth^2, data = .))$r.squared)) %>%
  ungroup()

p<-
ggplot() +
  geom_smooth(data = cog_df1, aes(x = cog_depth, y = temp), formula ='y ~ x',method = "lm", se = FALSE, color = 'black',alpha=0.7,linetype='dashed')+
  geom_point(data = na.omit(cog_df1), aes(x = cog_depth, y = temp, color = typeyear),alpha=0.6) +
  facet_wrap(region ~ sp, scales = 'free', nrow = 5) +
  labs(x = 'bathymetric center of gravity', y = 'mean region scaled SBT') +
  #scale_color_gradient2(low = 'blue',mid = 'white',high = 'red',midpoint = 0)+
  scale_color_manual(values = c('normal'='grey','cold'='blue',warm='red'),
                     breaks = unique(cog_df1[,c("typeyear","year")])[,'year'],  # Specify the breaks based on unique year values
                     labels = unique(cog_df1[,c("typeyear","year")])[,'year'],  # Use unique year values as labels
                     name = 'year') +
  theme_bw()+
  theme(aspect.ratio = 1,panel.grid.major = element_line(color = rgb(0, 0, 0,20, maxColorValue = 285), linetype = 'dashed', linewidth =  0.5),
        panel.background = element_rect(fill = NA),panel.ontop = TRUE,text = element_text(size=12),
        legend.background =  element_rect(fill = "transparent", colour = "transparent"),legend.key.height= unit(25, 'points'),
        legend.key.width= unit(25, 'points'),axis.title = element_blank(),legend.position = 'none',
        panel.border = element_rect(fill = NA, colour = 'black'),legend.key = element_rect(color="black"),
        axis.text = element_text(color='black'))+
  geom_text(data = rsquared_data, aes(label = sprintf("R^2 = %.2f", rsquared)), x = Inf, y = Inf, hjust = 1.1, vjust = 1.2)

  #geom_line(data = na.omit(cog_df1), aes(x = cog_depth, y = temp),stat="smooth",method = "lm", 
  #          linetype ="dashed",
  #          color='black')+

#save  plot
agg_png(paste0('./figures slope/regionSBT_cogdepth_full.png'), width = 20, height =13, units = "in", res = 300)
print(p)
dev.off()

#convert lat and lon for plotting reasons
cog_df2<-na.omit(cog_df1)
coordinates(cog_df2) <- ~ cog_lon + cog_lat
proj4string(cog_df2)<-CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
#proj4string(strata_pol) <- CRS('+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs') 
cog_df2<-spTransform(cog_df2,CRSobj = CRS('+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs'))
cog_df2<-as.data.frame(cog_df2,xy=TRUE)
cog_df2$year<-as.numeric(cog_df2$year)


#plot across simulations
p<-
ggplot()+
  geom_point(data=subset(cog_df2,typeyear=='normal'),aes(x=cog_lon,y=cog_lat,color=typeyear),size=0.7,alpha=0.9,shape=1)+
  geom_point(data=subset(cog_df2,typeyear %in% c('warm','cold')),aes(x=cog_lon,y=cog_lat,color=typeyear),size=0.7,alpha=0.7,shape=1)+
  geom_polygon(data=bs_sh,aes(x=long,y=lat,group=group),fill=NA,col='black')+
  geom_polygon(data=ak_sppoly,aes(x=long,y=lat,group=group),color='black',linewidth=0.2,fill = 'grey80')+
  #geom_point(data=cog_df1,aes(x=lon,y=lat,color=as.factor(year)))+
  #scale_color_manual(values = sunrise_palette, name='year')+
  #geom_path(data=cog_df2,aes(x=lon,y=lat))+
  scale_color_manual(values = c('normal'='grey','cold'='blue',warm='red'),
                     breaks = unique(cog_df2[,c("typeyear","year")])[,'year'],  # Specify the breaks based on unique year values
                     labels = unique(cog_df2[,c("typeyear","year")])[,'year'],  # Use unique year values as labels
                     name = 'year') +
  theme_bw()+
  theme(aspect.ratio = 1,panel.grid.major = element_line(color = rgb(0, 0, 0,20, maxColorValue = 285), linetype = 'dashed', linewidth =  0.5),
        panel.background = element_rect(fill = NA),panel.ontop = TRUE,text = element_text(size=12),
        legend.background =  element_rect(fill = "transparent", colour = "transparent"),legend.key.height= unit(25, 'points'),
        legend.key.width= unit(25, 'points'),axis.title = element_blank(),legend.position = c(0.12,0.47), #c(0.12,0.47)
        panel.border = element_rect(fill = NA, colour = 'black'),legend.key = element_rect(color="black"),
        axis.text = element_text(color='black'),legend.spacing.y = unit(10, 'points'),
        axis.text.y.right = element_text(hjust= 0.1 ,margin = margin(0,7,0,-25, unit = 'points'),color='black'),
        axis.text.x = element_text(vjust = 6, margin = margin(-7,0,7,0, unit = 'points'),color='black'),
        axis.ticks.length = unit(-5,"points"))+
  scale_x_continuous(expand = c(0,0),position = 'bottom',
                     breaks = c(-180,-170,-160),sec.axis = dup_axis())+
  scale_y_continuous(breaks = c(66,60,54))+
  coord_sf(crs = '+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs',
           xlim = panel_extent$x,
           ylim = c(panel_extent$y[1]-40000,panel_extent$y[2]-400000),
           label_axes = "-NE-")+
  labs(title=('CENTER OF GRAVITY by year'),xlab='',ylab='')+
  
  facet_wrap(region ~ sp, nrow = 5)
  #facet_grid(region ~ sp, scales = 'free', space = 'free', switch = 'y')

#save  plot
agg_png(paste0('./figures slope/regionSBT_coglonlat_full.png'), width = 20, height =13, units = "in", res = 300)
print(p)
dev.off()

#mean over simulations
cog_df3<-aggregate(cbind(cog_depth,cog_lon,cog_lat) ~ year+region+sp+temp+scaledtemp+typeyear,data=cog_df2,FUN=mean)

#plot mean over simulations
p<-
ggplot()+
  geom_point(data=subset(cog_df3,typeyear=='normal'),aes(x=cog_lon,y=cog_lat,color=typeyear),size=0.7,shape=19)+
  geom_point(data=subset(cog_df3,typeyear %in% c('warm','cold')),aes(x=cog_lon,y=cog_lat,color=typeyear),size=0.7,shape=19)+
  geom_polygon(data=bs_sh,aes(x=long,y=lat,group=group),fill=NA,col='black')+
  geom_polygon(data=ak_sppoly,aes(x=long,y=lat,group=group),color='transparent',linewidth=0.2,fill = 'grey80')+
  #geom_point(data=cog_df1,aes(x=lon,y=lat,color=as.factor(year)))+
  #scale_color_manual(values = sunrise_palette, name='year')+
  #geom_path(data=cog_df2,aes(x=lon,y=lat))+
  scale_color_manual(values = c('normal'='grey','cold'='blue',warm='red'),
                     breaks = unique(cog_df2[,c("typeyear","year")])[,'year'],  # Specify the breaks based on unique year values
                     labels = unique(cog_df2[,c("typeyear","year")])[,'year'],  # Use unique year values as labels
                     name = 'year') +
  theme_bw()+
  theme(aspect.ratio = 1,panel.grid.major = element_line(color = rgb(0, 0, 0,20, maxColorValue = 285), linetype = 'dashed', linewidth =  0.5),
        panel.background = element_rect(fill = NA),panel.ontop = TRUE,text = element_text(size=12),
        legend.background =  element_rect(fill = "transparent", colour = "transparent"),legend.key.height= unit(25, 'points'),
        legend.key.width= unit(25, 'points'),axis.title = element_blank(),legend.position = c(0.12,0.47), #c(0.12,0.47)
        panel.border = element_rect(fill = NA, colour = 'black'),legend.key = element_rect(color="black"),
        axis.text = element_text(color='black'),legend.spacing.y = unit(10, 'points'),
        axis.text.y.right = element_text(hjust= 0.1 ,margin = margin(0,7,0,-25, unit = 'points'),color='black'),
        axis.text.x = element_text(vjust = 6, margin = margin(-7,0,7,0, unit = 'points'),color='black'),
        axis.ticks.length = unit(-5,"points"))+
  scale_x_continuous(expand = c(0,0),position = 'bottom',
                     breaks = c(-180,-170,-160),sec.axis = dup_axis())+
  scale_y_continuous(breaks = c(66,60,54))+
  coord_sf(crs = '+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs',
           xlim = panel_extent$x,
           ylim = c(panel_extent$y[1]-40000,panel_extent$y[2]-400000),
           label_axes = "-NE-")+
  labs(title=('CENTER OF GRAVITY by year'),xlab='',ylab='')+
  facet_wrap(region ~ sp, nrow = 5)
#facet_grid(region ~ sp, scales = 'free', space = 'free', switch = 'y')

#save  plot
agg_png(paste0('./figures slope/regionSBT_coglonlat_mean.png'), width = 20, height =13, units = "in", res = 300)
print(p)
dev.off()

########################
########################
# correlation temp - depth
########################
########################

# # Pre-allocate memory for cog_df
# num_rows <- length(sp_shelfslope)*length(sims)*
#   ((length(dimnames(sim_dens1_slope)[[3]])*3)+(length(setdiff(yrs,dimnames(sim_dens1_slope)[[3]]))*5))
# 
# cog_df <- data.frame('sp' = character(num_rows),
#                      'sim' = character(num_rows),
#                      'year' = character(num_rows),
#                      'cog_lat' = numeric(num_rows),
#                      'cog_lon' = numeric(num_rows),
#                      'cog_depth' = numeric(num_rows),
#                      'region' = character(num_rows))
# 
# for (s in sp_sh) {
#   
# }
# 
# 
# # Calculate the correlation coefficient
# correlation <- cor(cog_df2$Temperature, data$Depth)
# 
# # Print the correlation coefficient
# print(correlation)
# 
# # Perform a correlation test
# cor.test(data$Temperature, data$Depth)


########################
########################
# Abundance fraction among regions
########################
########################

#df to store results
ind_df<-data.frame(matrix(NA,nrow = 0,ncol=4))
names(ind_df)<-c('sp','year','abundance','region')

#loop over common spp in shelf and slope
for (s in sp_shelfslope) {
  
  #s<-sp_shelfslope[1];y<-colnames(dens_slope)[1]
  
  cat(paste('#################### ',s," - ",'\n'))
  
  #EBS+NBS model fit and density
  load(paste0('./shelf EBS NBS VAST/',s,'/fit.RData'))
  fit_shelf<-fit
  bio_shelf<-fit_shelf$Report$Index_gctl[,,,1]
  
  #slope model fit ande density
  load(paste0('./slope EBS VAST/',s,'/fit.RData'))
  fit_slope<-fit
  bio_slope<-fit_slope$Report$Index_gctl[,,,1]
  
  #loop over years
  for (y in colnames(bio_slope)) {
    
    all_bio2<-c(bio_shelf[,y],bio_slope[,y])
    
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
      
      #append
      ind_df<-rbind(ind_df,
                    data.frame('sp'=s,
                               'year'=y,
                               'abundance'=sum(all_bio3),
                               'region'=r))  
      
    }
  }
}

# #save results
save(ind_df,file='./output/abundance.RData')
load(file='./output/abundance.RData')

#################

#filter for years in dimnames(sim_dens1_slope)[[3]]

#remove 'all' and 'EBS+NBS' region
ind_df1<-subset(ind_df,region %in% c('NBS','EBSshelf','EBSslope'))

ind_df2<-
ind_df1 %>% group_by(sp, year) %>%
  mutate(frac = abundance / sum(abundance))

ind_df2$year<-as.factor(ind_df2$year)

# Your original plot code with correction for coloring axis text
p<-
ggplot() +
  geom_bar(data = ind_df2, aes(x = year, y = frac, fill = region), stat = 'identity') +
  scale_fill_manual(values = c('EBSslope' = '#B4AF46', 'EBSshelf' = '#B4464B', 'NBS' = '#4682B4'),
                    breaks = c('EBSshelf', 'NBS', 'EBSslope'),
                    labels = c('EBSshelf', 'NBS', 'EBSslope')) +
  #scale_x_discrete(labels = ind_df2$year) +
  labs(y='abundance proportion')+
  theme_bw() +
  theme(text=element_text(size=14))+
  facet_wrap(~sp,scales='free_x', nrow = 5) +
  theme(axis.text.x = element_text(color = ifelse(levels(ind_df2$year) %in% c(2003:2005,2015:2016), "red", 
                                                  ifelse(levels(ind_df2$year) %in% c(2009,2008,2012), "blue",'black'))))

#save  plot
agg_png(paste0('./figures slope/abundance_proportion.png'), width = 14, height =10, units = "in", res = 300)
print(p)
dev.off()


########################
########################
########################
########################

#adding outer strata
#check and compare fraction and cog depth midshelf outershelf

#######################
#######################

#maps add coldpool

#add mean depth 

#only for years with slope survey

#slow abundance years vs high abundance years (bubble plot to add another variable but we could remove map)

#summary of species responses - expand?
#pollock an cod respond and use the space different from arrowtooth flounder

#for slope we could predict over different depth by year 

#compare outer shelf and slope?

#calculate edges 95% abundance by depth relative to the mean


#######################
########################
# Abundance fraction among depth regions
########################
########################

#df to store results
ind_df<-data.frame(matrix(NA,nrow = 0,ncol=5))
names(ind_df)<-c('sp','year','abundance','dens','region')

#loop over common spp in shelf and slope
for (s in sp_shelfslope) {
  
  #s<-sp_shelfslope[1];y<-colnames(dens_slope)[1]
  
  cat(paste('#################### ',s," - ",'\n'))
  
  #EBS+NBS model fit and density
  load(paste0('./shelf EBS NBS VAST/',s,'/fit.RData'))
  fit_shelf<-fit
  bio_shelf<-fit_shelf$Report$Index_gctl[,,,1]
  
  #slope model fit ande density
  load(paste0('./slope EBS VAST/',s,'/fit.RData'))
  fit_slope<-fit
  bio_slope<-fit_slope$Report$Index_gctl[,,,1]
  
  #loop over years
  for (y in colnames(bio_slope)) {
    
    all_bio2<-c(bio_shelf[,y],bio_slope[,y])
    
    #loop over regions
    for (r in c(na.omit(unique(grid.ebs_year1$region_depth)),'200_400','0_200')) {
      
      if (r=='0_100') {
        c<-which(grid.ebs_year1$region_depth=="0_100")
      } else if (r=='100_200') {
        c<-which(grid.ebs_year1$region_depth=='100_200')
      } else if (r == '200_300') {
        c <- which(grid.ebs_year1$region_depth=='200_300')
      } else if (r=="300_400") {
        c<-which(grid.ebs_year1$region_depth=='300_400')
      } else if (r=='200_400') {
        c<-which(grid.ebs_year1$region_depth %in% c('200_300','300_400'))
      } else if (r=='0_200') {
        c<-which(grid.ebs_year1$region_depth %in% c('0_100','100_200'))
      }
      
      #subset grid and dens by grid
      igrid<-grid.ebs_year1[c,]
      all_bio3<-all_bio2[c]
      
      #calculate percentiles
      #data.frame(t(quantile(all_bio3, probs = c(0.05, 0.10, 0.90, 0.95)))) # Calculate percentiles
      
      #bio to dens
      all_dens3<-all_bio2[c]/igrid$Area_in_survey_km2      
      
      #append
      ind_df<-rbind(ind_df,
                    data.frame('sp'=s,
                               'year'=y,
                               'abundance'=sum(all_bio3),
                               'dens'=mean(all_dens3),
                               'region'=r,
                               'percentiles'=data.frame(t(quantile(all_bio3, probs = c(0.05, 0.10, 0.90, 0.95)))))) # Calculate percentiles)))  
      
      
      
    }
  }
}

# #save results
save(ind_df,file='./output/abundance_depth.RData')
load(file='./output/abundance_depth.RData')

#subset
#remove 'all' and 'EBS+NBS' region
#ind_df1<-subset(ind_df,region %in% c('0_100','100_200','200_300','300_400'))
ind_df1<-subset(ind_df,region %in% c('0_200','200_400'))
ind_df1<-subset(ind_df,region %in% c('0_100','100_200'))

##########################
# RANGE BY DEPTH AREAS
#########################

ind_df2<-
  ind_df1 %>% group_by(sp, year) %>%
  mutate(frac = abundance / sum(abundance))

ind_df2$year<-as.factor(ind_df2$year)


#plot
#p<-
ggplot() +
  geom_errorbar(data = ind_df1, aes(x = year, ymax=percentiles.X95.,ymin=percentiles.X5., color = region),position = position_dodge(width = .5)) +
  # scale_fill_manual(values = c('EBSslope' = '#B4AF46', 'EBSshelf' = '#B4464B', 'NBS' = '#4682B4'),
  #                   breaks = c('EBSshelf', 'NBS', 'EBSslope'),
  #                   labels = c('EBSshelf', 'NBS', 'EBSslope')) +
  #scale_x_discrete(labels = ind_df2$year) +
  labs(y='range biomass (5% - 95%)')+
  theme_bw() +
  theme(text=element_text(size=14))+
  facet_wrap(~sp,scales='free', nrow = 5) +
  theme(axis.text.x = element_text(color = ifelse(levels(ind_df2$year) %in% c(2003:2005,2015:2016), "red", 
                                                  ifelse(levels(ind_df2$year) %in% c(2009,2008,2012), "blue",'black'))))



##########################
# ABUNDANCE PROPORTION
#########################

#plot
#p<-
  ggplot() +
  geom_bar(data = ind_df2, aes(x = year, y = frac, fill = region), stat = 'identity') +
  # scale_fill_manual(values = c('EBSslope' = '#B4AF46', 'EBSshelf' = '#B4464B', 'NBS' = '#4682B4'),
  #                   breaks = c('EBSshelf', 'NBS', 'EBSslope'),
  #                   labels = c('EBSshelf', 'NBS', 'EBSslope')) +
  #scale_x_discrete(labels = ind_df2$year) +
  labs(y='abundance proportion')+
  theme_bw() +
  theme(text=element_text(size=14))+
  facet_wrap(~sp,scales='free_x', nrow = 5) +
  theme(axis.text.x = element_text(color = ifelse(levels(ind_df2$year) %in% c(2003:2005,2015:2016), "red", 
                                                  ifelse(levels(ind_df2$year) %in% c(2009,2008,2012), "blue",'black'))))
  
  
##########################
# MEAN DENSITY
######################### 

#plot
p<-
  ggplot() +
    geom_point(data = ind_df1, aes(x = year, y = dens, color = region), stat = 'identity') +
    # scale_fill_manual(values = c('EBSslope' = '#B4AF46', 'EBSshelf' = '#B4464B', 'NBS' = '#4682B4'),
    #                   breaks = c('EBSshelf', 'NBS', 'EBSslope'),
    #                   labels = c('EBSshelf', 'NBS', 'EBSslope')) +
    #scale_x_discrete(labels = ind_df2$year) +
    labs(y='mean density')+
    theme_bw() +
    theme(text=element_text(size=14))+
    facet_wrap(~sp,scales='free_y', nrow = 5) +
    theme(axis.text.x = element_text(color = ifelse(levels(ind_df2$year) %in% c(2003:2005,2015:2016), "red", 
                                                    ifelse(levels(ind_df2$year) %in% c(2009,2008,2012), "blue",'black'))))
  

#save  plot
agg_png(paste0('./figures slope/abundance_proportion.png'), width = 14, height =10, units = "in", res = 300)
print(p)
dev.off()

##########################
# RANGE EDGES BY DEPTH AREAS
#########################

#grid without slope
grid.ebs_year1shelf<-subset(grid.ebs_year1,region %in% c('EBSshelf','NBS','EBSslope'))

#loop over common spp in shelf and slope
for (s in sp_shelfslope) {
  
  s<-sp_shelfslope[1]#;y<-colnames(dens_slope)[1]
  
  
  #EBS+NBS model fit and density
  load(paste0('./shelf EBS NBS VAST/',s,'/fit.RData'))
  fit_shelf<-fit
  bio_shelf<-fit_shelf$Report$Index_gctl[,,,1]
  
  #slope model fit ande density
  load(paste0('./slope EBS VAST/',s,'/fit_st.RData'))
  fit_slope<-fit
  bio_slope<-fit_slope$Report$Index_gctl[,,,1] #biomass
  
  plot_list<-list()
  
  #loop over years
  for (y in as.factor(2002:2016)) {
    
    #y<-'2016'
    
    cat(paste('#################### ',s," - ",y,'\n'))
    
    
    #grid without slope
    grid.ebs_year1<-subset(grid.ebs_year,Year==y)
    
    all_bio2<-c(bio_shelf[,y],bio_slope[,y])
    all_bio3<-cbind(all_bio2,grid.ebs_year1)
    all_bio3$dens<-all_bio3$all_bio2/all_bio3$Area_in_survey_km2
    all_bio4<-subset(all_bio3,region %in% c('EBSshelf','EBSslope','NBS') & DepthGEBCO >0)
    
    p5<-quantile(all_bio4$dens,c(0.10))
    p95<-quantile(all_bio4$dens,c(0.90))
    
    
    all_bio4$quantile<-ifelse(all_bio4$dens < p95 & all_bio4$dens > p5,1,0)
    all_bio4$all<-1
    
    #df to spatialpoint df
    coordinates(all_bio4) <- ~ Lon + Lat
    crs(all_bio4)<-c(crs="+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
    
    #reproject coordinates for plotting purposes
    all_bio5<-spTransform(all_bio4,'+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs')
    all_bio6<-data.frame(all_bio5)
    
    #x and y cells
    xycells<-as.integer(sqrt(dim(all_bio5)[1]))
    
    # create a template raster
    r1 <- raster(ext=extent(all_bio5),ncol=xycells, nrow=xycells) #c(15800,15800) 7000
    
    #create raster
    r2<-rasterize(all_bio5, r1 ,field=c('quantile','Temp','all'))
    crs(r2) <- CRS('+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs')
    r2$quantile[r2$quantile==0]<-NA
    r2$Temp[r2$Temp>2]<-NA
    r2[r2$Temp<2]<-1
    
    
    #create polygon to get boundaries of each strata
    r3<-as.data.frame(r2,xy=TRUE)
    r3<-r3[complete.cases(r3$quantile),] 
    #as factors
    r3$quantile<-as.factor(r3$quantile)
    r4a<-rasterToPolygons(r2$quantile,dissolve=TRUE,digits = 1)
    r4b<-rasterToPolygons(r2$Temp,dissolve=TRUE,digits = 1)
    r4c<-rasterToPolygons(r2$all,dissolve=TRUE,digits = 1)
    
    # Function to get color based on temperature
    get_temp_color <- function(itemp, temp_min, temp_max) {
      gradient <- colorRampPalette(c("blue", "white", "darkred"))(100)
      temp_scaled <- scales::rescale(itemp, to = c(1, 100), from = c(temp_min, temp_max))
      return(gradient[round(temp_scaled)])
    }
    
    # Define the temperature range
    alltemp_region<-subset(temp_region,region=='all')
    temp_min <- min(alltemp_region$temp)
    temp_max <- max(alltemp_region$temp)
    # Extract temperature for the year y
    itemp_df <- subset(alltemp_region, year == y)
    itemp <- itemp_df$temp
    
    #plot
    p<-
    ggplot()+
      geom_polygon(data=r4b,aes(x=long,y=lat,group=group,color='coldpool'), fill='blue',alpha=0.7)+
      geom_tile(data=r3,aes(x=x,y=y,fill=quantile),alpha=0.4)+
      scale_fill_manual(values='grey30',labels='species range')+
      geom_polygon(data=r4a,aes(x=long,y=lat,group=group),fill='transparent', colour="grey50")+
      scale_color_manual(values=c('cold pool'='transparent'))+
      geom_polygon(data=ak_sppoly,aes(x=long,y=lat,group=group),color='transparent',linewidth=0.2,fill = 'white')+
      geom_polygon(data=NBS_sh,aes(x=long,y=lat,group=group),fill=NA,col='black')+
      geom_polygon(data=EBSshelf_sh,aes(x=long,y=lat,group=group),fill=NA,col='black')+
      geom_polygon(data=EBSslope_sh,aes(x=long,y=lat,group=group),fill=NA,col='black')+
      theme_bw()+
      labs(title = y)+
      annotate("text", x = -1200000, y = 1750000, label = y, size = 10, color = get_temp_color(itemp, temp_min, temp_max)) +
      coord_sf(crs = '+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs',
              #xlim = c(panel_extent$x[1]+200000,panel_extent$x[2]),
              xlim = c(panel_extent$x[1]+200000,panel_extent$x[2]+30000),
              ylim = c(panel_extent$y[1]-100000,panel_extent$y[2]-200000),
              label_axes = "-NE-")+
      theme(axis.title = element_blank(),
            legend.title = element_blank(),
            plot.margin = unit(c(0, 0, 0, 0),units = 'in'), 
            axis.text = element_blank(),
            legend.position = 'none',
            panel.grid = element_blank(),
            plot.title = element_blank(),
            axis.ticks = element_blank())
    
    plot_list[[paste0(s, "_", y)]] <- p
    
    # #loop over regions
    # for (r in c(na.omit(unique(grid.ebs_year1$region_depth)),'200_400','0_200')) {
    #   
    #   if (r=='0_100') {
    #     c<-which(grid.ebs_year1$region_depth=="0_100")
    #   } else if (r=='100_200') {
    #     c<-which(grid.ebs_year1$region_depth=='100_200')
    #   } else if (r == '200_300') {
    #     c <- which(grid.ebs_year1$region_depth=='200_300')
    #   } else if (r=="300_400") {
    #     c<-which(grid.ebs_year1$region_depth=='300_400')
    #   } else if (r=='200_400') {
    #     c<-which(grid.ebs_year1$region_depth %in% c('200_300','300_400'))
    #   } else if (r=='0_200') {
    #     c<-which(grid.ebs_year1$region_depth %in% c('0_100','100_200'))
    #   }
    #   
    #   #subset grid and dens by grid
    #   igrid<-grid.ebs_year1[c,]
    #   all_bio3<-all_bio2[c]
    #   
    #   #calculate percentiles
    #   #data.frame(t(quantile(all_bio3, probs = c(0.05, 0.10, 0.90, 0.95)))) # Calculate percentiles
    #   
    #   #bio to dens
    #   all_dens3<-all_bio2[c]/igrid$Area_in_survey_km2      
    #   
    #   #append
    #   ind_df<-rbind(ind_df,
    #                 data.frame('sp'=s,
    #                            'year'=y,
    #                            'abundance'=sum(all_bio3),
    #                            'dens'=mean(all_dens3),
    #                            'region'=r,
    #                            'percentiles'=data.frame(t(quantile(all_bio3, probs = c(0.05, 0.10, 0.90, 0.95)))))) # Calculate percentiles)))  
      
      
  }
  
 
  # Combine all plots into a single data frame
  plot_data <- do.call(rbind, lapply(names(plot_list), function(name) {
    p <- plot_list[[name]]
    data.frame(year = substr(name, nchar(s) + 2, nchar(name)), plot = I(list(p)))
  }))
  
  # Create the animation
  animated_plot <- ggplot(plot_data, aes(x = x, y = y)) +
    geom_tile(data = plot_data, aes(fill = plot), alpha = 0.6) +
    transition_states(year, transition_length = 2, state_length = 1) +
    ease_aes('linear') +
    labs(title = "Year: {closest_state}") +
    theme_minimal()
  
  
  
  
  p<-plot_list[['Gadus chalcogrammus_2002']]
  year <- substr('Gadus chalcogrammus_2002', nchar(s) + 2, nchar('Gadus chalcogrammus_2002'))
  data.frame(year = year, plot = I(list(p)))
  
  
  legend1 <- cowplot::get_legend( 
    plot_list[[1]] + 
      theme(legend.position = "right",legend.spacing.y = unit(-0.1,'in')) #,legend.justification = c(0.72,0.5)
  ) 
  
  #plot(legend1)
  
  
  plot_list[['2023']]<-legend1
  
  agg_png(paste0('./figures slope/',s,'_rangeEBS.png'), width = 10, height = 10, units = "in", res = 300)
  print(
  plot_grid(plotlist = plot_list)
  )
  dev.off()
  
}

##########################
# USING EBS DEPTH RANGES as 0-50 coastal 50-100 middle 100-200 outer
#########################

#add depths regions
region_depth1<-ifelse(grid.ebs_year1$DepthGEBCO<50,'coastal',
                     ifelse(grid.ebs_year1$DepthGEBCO>=50 & grid.ebs_year1$DepthGEBCO<=100,'middle',
                            ifelse(grid.ebs_year1$DepthGEBCO>=100 & grid.ebs_year1$DepthGEBCO<=200,'outer',NA)))

region_depth1<-ifelse(grid.ebs_year1$region=='NBS',NA,region_depth1)

#cbind
grid.ebs_year1<-cbind(grid.ebs_year1,region_depth1)
grid.ebs_year2<-subset(grid.ebs_year1,region=='EBSshelf')

#df to store results
ind_df1<-data.frame(matrix(NA,nrow = 0,ncol=5))
names(ind_df1)<-c('sp','year','abundance','dens','region')

#loop over common spp in shelf and slope
for (s in sp_shelfslope) {
  
  #s<-sp_shelfslope[1];y<-colnames(dens_slope)[1]
  
  cat(paste('#################### ',s," - ",'\n'))
  
  #EBS+NBS model fit and density
  load(paste0('./shelf EBS NBS VAST/',s,'/fit.RData'))
  fit_shelf<-fit
  bio_shelf<-fit_shelf$Report$Index_gctl[,,,1]
  
  #slope model fit ande density
  load(paste0('./slope EBS VAST/',s,'/fit.RData'))
  fit_slope<-fit
  bio_slope<-fit_slope$Report$Index_gctl[,,,1]
  
  #loop over years
  for (y in colnames(bio_shelf)) {
    
    all_bio2<-c(bio_shelf[,y])
    
    #loop over regions
    for (r in c(na.omit(unique(grid.ebs_year1$region_depth1)))) {
      
      if (r=='coastal') {
        c<-which(grid.ebs_year2$region_depth1=="coastal")
      } else if (r=='middle') {
        c<-which(grid.ebs_year2$region_depth1=='middle')
      } else if (r == 'outer') {
        c <- which(grid.ebs_year2$region_depth1=='outer')
      } 
      
      #subset grid and dens by grid
      igrid<-grid.ebs_year2[c,]
      all_bio3<-all_bio2[c]
      
      #calculate percentiles
      #data.frame(t(quantile(all_bio3, probs = c(0.05, 0.10, 0.90, 0.95)))) # Calculate percentiles
      
      #bio to dens
      all_dens3<-all_bio2[c]/igrid$Area_in_survey_km2      
      
      #append
      ind_df1<-rbind(ind_df1,
                    data.frame('sp'=s,
                               'year'=y,
                               'abundance'=sum(all_bio3),
                               'dens'=mean(all_dens3),
                               'region'=r,
                               'percentiles'=data.frame(t(quantile(all_bio3, probs = c(0.05, 0.10,0.20,0.80, 0.90, 0.95)))))) # Calculate percentiles)))  
      
      
      
    }
  }
}

# #save results
save(ind_df1,file='./output/abundance_depth1.RData')
load(file='./output/abundance_depth1.RData')


##########################
# RANGE BY DEPTH AREAS
#########################

ind_df2<-
  ind_df1 %>% group_by(sp, year) %>%
  mutate(frac = abundance / sum(abundance))

ind_df2$year<-as.factor(ind_df2$year)

ind_df3<-subset(ind_df2,year %in% c(2016,2018,2019,2022,1992,2019,2007:2010,2012))
ind_df3$year<-as.factor(ind_df3$year)
# Drop unused levels
ind_df3$year <- droplevels(ind_df3$year)

#plot
#p<-
ggplot() +
  geom_errorbar(data = ind_df3, aes(x = year, ymax=percentiles.X95.,ymin=percentiles.X5., color = region),position = position_dodge(width = .5)) +
  # scale_fill_manual(values = c('EBSslope' = '#B4AF46', 'EBSshelf' = '#B4464B', 'NBS' = '#4682B4'),
  #                   breaks = c('EBSshelf', 'NBS', 'EBSslope'),
  #                   labels = c('EBSshelf', 'NBS', 'EBSslope')) +
  #scale_x_discrete(labels = ind_df2$year) +
  labs(y='range biomass (10% - 90%)')+
  #scale_x_discrete(breaks=c(2016,2018,2019,2022,1992,2019,2007:2010,2012))+
  theme_bw() +
  theme(text=element_text(size=14))+
  facet_wrap(~sp,scales='free', nrow = 5) +
  theme(axis.text.x = element_text(color = ifelse(levels(ind_df3$year) %in% c(2016,2018,2019,2022), "red",
                                                  ifelse(levels(ind_df3$year) %in% c(1992,2007:2010,2012), "blue",'black'))))
  # theme(axis.text.x = element_text(color = ifelse(levels(ind_df2$year) %in% c(2003:2005,2015:2016), "red",
  #                                                 ifelse(levels(ind_df2$year) %in% c(2009,2008,2012), "blue",'black'))))



##########################
# ABUNDANCE PROPORTION
#########################

#plot
#p<-
ggplot() +
  geom_bar(data = ind_df3, aes(x = year, y = frac, fill = region), stat = 'identity') +
  # scale_fill_manual(values = c('EBSslope' = '#B4AF46', 'EBSshelf' = '#B4464B', 'NBS' = '#4682B4'),
  #                   breaks = c('EBSshelf', 'NBS', 'EBSslope'),
  #                   labels = c('EBSshelf', 'NBS', 'EBSslope')) +
  #scale_x_discrete(labels = ind_df2$year) +
  labs(y='abundance proportion')+
  theme_bw() +
  theme(text=element_text(size=14))+
  facet_wrap(~sp,scales='free_x', nrow = 5) +
  theme(axis.text.x = element_text(color = ifelse(levels(ind_df3$year) %in% c(2016,2018,2019,2022), "red",
                                                  ifelse(levels(ind_df3$year) %in% c(1992,2007:2010,2012), "blue",'black'))))
# theme(axis.text.x = element_text(color = ifelse(levels(ind_df2$year) %in% c(2003:2005,2015:2016), "red",
#                                                 ifelse(levels(ind_df2$year) %in% c(2009,2008,2012), "blue",'black'))))


##########################
# MEAN DENSITY
#########################
ind_df11<-
  
  # Scaling values by group
  ind_df11 <- data.frame(ind_df1 %>%
                               group_by(sp,region) %>%
                               mutate(scaleddens = scale_value(dens)) %>%
                               ungroup())
  

#plot
p<-
  ggplot() +
  geom_point(data = ind_df11, aes(x = year, y = scaleddens, color = region), stat = 'identity') +
  # scale_fill_manual(values = c('EBSslope' = '#B4AF46', 'EBSshelf' = '#B4464B', 'NBS' = '#4682B4'),
  #                   breaks = c('EBSshelf', 'NBS', 'EBSslope'),
  #                   labels = c('EBSshelf', 'NBS', 'EBSslope')) +
  #scale_x_discrete(labels = ind_df2$year) +
  labs(y='mean density')+
  theme_bw() +
  theme(text=element_text(size=14))+
  facet_wrap(~sp,scales='free_y', nrow = 5) +
  theme(axis.text.x = element_text(color = ifelse(levels(ind_df2$year) %in% c(2003:2005,2015:2016), "red", 
                                                  ifelse(levels(ind_df2$year) %in% c(2009,2008,2012), "blue",'black'))))


#save  plot
agg_png(paste0('./figures slope/abundance_proportion.png'), width = 14, height =10, units = "in", res = 300)
print(p)
dev.off()


########################################
#
# PREDICT SLOPE MODEL OVER RANGE OF DEPTH FOR EVERY YEAR
#
#######################################

library(VAST)
library(effects)  # Used to visualize covariate effects
library(splines)

scale_value <- function(x) {
  (x - mean(x)) / sd(x)
}


sp<-sp_shelfslope[1]

load(paste0('./slope EBS VAST/',sp,'/fit.RData'))

df1<-readRDS(paste0('./data processed/species/',sp,'/data_geostat.rds'))
df2<-subset(df1,survey_name==unique(df1$survey_name)[3])
df2$presence<-ifelse(df2$cpue_kgha!=0,1,0)
df2$ScaleLogDepth<-scale_value(df2$depth_m)

df2$year_type<-ifelse(df2$year %in% c(2016,2018,2019,2022), "warm",
                      ifelse(df2$year %in% c(2002,2004),'normal','cold'))

year_colors <- c("2002" = "black", "2004" = "black", "2008" = "blue", "2010" = "blue", "2012" = "blue", "2016" = "red")


ggplot()+
  geom_point(data =subset(df2,cpue_kgha!=0), aes(x = ScaleLogDepth, y = log(cpue_kgha),color=as.factor(year)),alpha=0.6) +
  labs(x = 'scaled log depth', y = 'log positive biomass',title = sp) +
  theme_bw()+
  scale_color_manual(values = year_colors) +
  #theme(legend.position = 'none')+
  geom_smooth(data = subset(df2,cpue_kgha!=0), aes(x = ScaleLogDepth, y = log(cpue_kgha),color=as.factor(year),group=as.factor(year)), method = "lm", se = FALSE,alpha=0.9, 
              formula= 'y ~ poly(x, 2)')


# geom_line(data = obs1, aes(x = dens, y = ScaleLogDepth),stat="smooth",method = "lm", 
#           linetype ="dashed",
#           color='black')+
#scale_size() #limits = c(100000,1000000000000),



#fit
# Must add data-frames to global environment (hope to fix in future)
covariate_data_full = fit$effects$covariate_data_full
catchability_data_full = fit$effects$catchability_data_full

# Plot 1st linear predictor, but could use `transformation` to apply link function
pred = Effect.fit_model( fit,
                         focal.predictors = c("ScaleLogDepth"),
                         which_formula = "X2",
                         xlevels = 100,
                         transformation = list(link=identity, inverse=identity) )

plot(pred)

ggplot()+
  geom_point(data =df2, aes(x = ScaleLogDepth, y = presence,color=as.factor(year)),alpha=0.6) +
  #facet_wrap( ~ year, scales = 'free', nrow = 5) +
  labs(x = 'depth', y = 'positive biomass') +
  theme_bw()+
  theme(legend.position = 'none')+
  geom_smooth(data = df2, aes(x = ScaleLogDepth, y = presence,color=as.factor(year),group=as.factor(year)), method = "loess", se = FALSE,alpha=0.9)

#fit
# Must add data-frames to global environment (hope to fix in future)
covariate_data_full = fit$effects$covariate_data_full
catchability_data_full = fit$effects$catchability_data_full

# Plot 1st linear predictor, but could use `transformation` to apply link function
pred = Effect.fit_model( fit,
                         focal.predictors = c("ScaleLogDepth"),
                         which_formula = "X1",
                         xlevels = 100,
                         transformation = list(link=identity, inverse=identity) )

plot(pred)




grid_slope$ScaleLogDepth<-scale_value(log(grid_slope$DepthGEBCO))
newdata<-grid_slope[,c("Year","Lat","Lon","ScaleLogDepth",'DepthGEBCO')]
names(newdata)[5]<-'Depth'
newdata$CPUE_kg<-0
#newdata$Year<-NA

# library(pdp)
# 
# # Make function to interface with pdp
# pred.fun = function( object, newdata ){
#   predict( x=object,
#            Lat_i = object$data_frame$Lat_i,
#            Lon_i = object$data_frame$Lon_i,
#            t_i = object$data_frame$t_i,
#            a_i = object$data_frame$a_i,
#            what = "P1_iz",
#            new_covariate_data = newdata,
#            do_checks = FALSE )
# }
# 
# names(fit$covariate_data);summary(fit$covariate_data)
# names(newdata);summary(newdata)
# # Run partial
# Partial = partial( object = fit,
#                    pred.var = "ScaleLogDepth",
#                    pred.fun = pred.fun,
#                    train = fit$covariate_data )


library(effects)  # Used to visualize covariate effects
library(splines)

#fit
# Must add data-frames to global environment (hope to fix in future)
covariate_data_full = fit$effects$covariate_data_full
catchability_data_full = fit$effects$catchability_data_full

# Plot 1st linear predictor, but could use `transformation` to apply link function
pred = Effect.fit_model( fit,
                         focal.predictors = c("ScaleLogDepth"),
                         which_formula = "X2",
                         xlevels = 100,
                         transformation = list(link=identity, inverse=identity) )

plot(pred)
#
obs<-
cbind(
fit$data_frame,
'ScaleLogDepth'=fit$covariate_data$ScaleLogDepth)
obs$dens<-obs$b_i/obs$a_i
obs1<-obs[pred_TF == 1,]


#df<-
  Sim1$b_i

ggplot()+
  geom_point(data = subset(obs1,b_i!=0), aes(x = ScaleLogDepth, y = dens),alpha=0.6) +
  facet_wrap( ~ t_i, scales = 'free', nrow = 5) +
  labs(x = 'bathymetric center of gravity', y = 'mean region SBT') +
  theme_bw()+
  theme(legend.position = 'none')#+
  #geom_smooth(data = cog_df1, aes(x = depth, y = temp), method = "lm", se = FALSE, color = 'black',alpha=0.9)
  # geom_line(data = obs1, aes(x = dens, y = ScaleLogDepth),stat="smooth",method = "lm", 
  #           linetype ="dashed",
  #           color='black')+
  #scale_size() #limits = c(100000,1000000000000),


#fit1
# Must add data-frames to global environment (hope to fix in future)
covariate_data_full = fit1$effects$covariate_data_full
catchability_data_full = fit1$effects$catchability_data_full

# Plot 1st linear predictor, but could use `transformation` to apply link function
pred = Effect.fit_model( fit1,
                         focal.predictors = c("ScaleLogDepth"),
                         which_formula = "X2",
                         xlevels = 100,
                         transformation = list(link=identity, inverse=identity) )



#predictions and data in histograms by bins of depth / slope vs outer (slope: 200-400; outer: 150-200) _ DONE

#change transparency of blue polygon _ DONE look the same

#change the percentile80-20%? and check upper seperately _DONE

#check max depth  _DONE

#mean cog depth warm-cold _DONE

#get temperature of EBS and wamr/cold years when evaluating slope

#cod and pollock should show that shifts distribution

#effective area occupied Thorson et al 2016
#ht = bt / mt
#area occupied = total abundance / average density


############################################################################################################################################################
############################################################################################################################################################
############################################################################################################################################################
############################################################################################################################################################
############################################################################################################################################################

################
# DATA
################

for (s in sp_shelfslope) {
  
  #s<-sp_shelfslope[2]
  
  #####################
  # DENSITY DATA
  #####################
  
  #shelf
  data<-readRDS(paste0('./data processed/species/',s,'/data_geostat_temp.rds'))
  
  #filter by region and years
  data_shelf1<-data[which(data$survey_name=='Eastern Bering Sea Crab/Groundfish Bottom Trawl Survey'),]
  data_shelf2<-data_shelf1[which(data_shelf1$year %in% c(1992,2007:2010,2012,2016,2018,2019,2022)),]
  data_shelf3<-data_shelf2[,c("lat_start","lon_start","year", "cpue_kgkm2" , "depth_m")]
  
  #slope
  data_slope1<-data[which(data$survey_name=='Eastern Bering Sea Slope Bottom Trawl Survey'),]
  data_slope2<-data_slope1[which(data_slope1$year %in% c(1992,2007:2010,2012,2016,2018,2019,2022)),]
  data_slope3<-data_slope2[,c("lat_start","lon_start","year", "cpue_kgkm2" , "depth_m")]
  #EBS_slope: 2002, 2004, 2008, 2010, 2012 and 2016 
  
  data<-rbind(data_shelf3,data_slope3)
  
  data$depth_bin<-ifelse(data$depth_m<50,'0-50',
                         ifelse(data$depth_m<100,'50-100',
                                ifelse(data$depth_m<150,'100-150',
                                       ifelse(data$depth_m<200,'150-200','>200'))))
                                       
  data$year_type<-ifelse(data$year %in% c(2016,2018,2019,2022), "warm",'cold')
  data$year_type<-as.factor(data$year_type)                
  
  data$depth_bin<-factor(data$depth_bin, levels = c("0-50", "50-100", '100-150',"150-200",">200"))
  data$year<-as.factor(data$year)
  
  #plot
  p1<-
  ggplot()+
    geom_boxplot(data=data,aes(x=year,y=log(cpue_kgkm2),fill=depth_bin,group=interaction(year,depth_bin)),position=position_dodge(width=0.8,preserve = "single"),outliers = FALSE)+
    scale_fill_brewer(name = 'depth bins (m)')+
    theme_bw()+
    theme(axis.title.x = element_blank(),axis.text.x = element_text(color = ifelse(levels(data$year) %in% c(2016,2018,2019,2022), "red",
                                                  ifelse(levels(data$year) %in% c(1992,2007:2010,2012), "blue",'black')),size=14),plot.title = element_text(size=16),legend.text = element_text(size = 12),  # Increase legend text size
          legend.title = element_text(size = 14), # Increase legend title size
          legend.key.size = unit(1.5, 'lines'),axis.title.y = element_text(size=14),
          legend.position = 'none')+
    ggtitle(label =s)+
    ylab(label = 'log(density) (kg/km2)')#+
  
  p2<-
    ggplot()+
      geom_boxplot(data=data,aes(x=year_type,y=log(cpue_kgkm2),fill=depth_bin,group=interaction(year_type,depth_bin)),position=position_dodge(width=0.8,preserve = "single"),outliers = FALSE)+
      scale_fill_brewer(name = 'depth bins (m)')+
      theme_bw()+
      theme(axis.title.x = element_blank(),axis.text.x = element_text(color = ifelse(levels(data$year_type) == 'warm', "red",'blue'),size=14),axis.title.y = element_blank(),plot.title = element_text(size=16),legend.text = element_text(size = 12),  # Increase legend text size
            legend.title = element_text(size = 14), # Increase legend title size
            legend.key.size = unit(1.5, 'lines'))+
      ggtitle(label ='')#+
      #ylab(label = 'log(density) (kg/km2)')#+
  
  # print(
  # cowplot::plot_grid(p1,p2,rel_widths = c(1,0.6))
  # )
  
  #####################
  # PERCENTILE DATA
  #####################
  
  data_nonzero<-subset(data,cpue_kgkm2!=0)
  
  per_data<-
  aggregate(cpue_kgkm2 ~ year + depth_bin,data_nonzero,quantile,probs=0.2)
  
  per_data$year_type<-ifelse(per_data$year %in% c(2016,2018,2019,2022), "warm",'cold')
  per_data$year_type<-as.factor(per_data$year_type)                
  
  per_data$depth_bin<-factor(per_data$depth_bin, levels = c("0-50", "50-100", '100-150',"150-200",">200"))
  per_data$year<-as.factor(per_data$year)
  
  p1<-
  ggplot()+
    geom_point(data=per_data,aes(x=year,y=log(cpue_kgkm2),color=depth_bin,group=interaction(year,depth_bin)),position=position_dodge(width=0.5))+
    scale_color_brewer(name = 'depth bins (m)')+
    theme_bw()+
    theme(axis.title.x = element_blank(),axis.text.x = element_text(color = ifelse(levels(data$year) %in% c(2016,2018,2019,2022), "red",
                                                                                   ifelse(levels(data$year) %in% c(1992,2007:2010,2012), "blue",'black')),size=14),plot.title = element_text(size=16),legend.text = element_text(size = 12),  # Increase legend text size
          legend.title = element_text(size = 14), # Increase legend title size
          legend.key.size = unit(1.5, 'lines'),axis.title.y = element_text(size=14),
          legend.position = 'none')+
    ggtitle(label =s)+
    ylab(label = 'percentile 20% log(density)')#+
  
  p2<-
    ggplot()+
    geom_boxplot(data=per_data,aes(x=year_type,y=log(cpue_kgkm2),fill=depth_bin,group=interaction(year_type,depth_bin)),position=position_dodge(width=0.8,preserve = "single"),outliers = FALSE)+
    scale_fill_brewer(name = 'depth bins (m)')+
    theme_bw()+
    theme(axis.title.x = element_blank(),axis.text.x = element_text(color = ifelse(levels(data$year_type) == 'warm', "red",'blue'),size=14),axis.title.y = element_blank(),plot.title = element_text(size=16),legend.text = element_text(size = 12),  # Increase legend text size
          legend.title = element_text(size = 14), # Increase legend title size
          legend.key.size = unit(1.5, 'lines'))+
    ggtitle(label ='')#+
    #ylab(label = 'log(density) (kg/km2)')#+
  
  print(
    cowplot::plot_grid(p1,p2,rel_widths = c(1,0.6))
  )

  #####################
  # COG DEPTH
  #####################
  
  data$dens_depth<-data$cpue_kgkm2*data$depth_m
  
  cog_data<-
  aggregate(cbind(dens_depth,cpue_kgkm2) ~ year ,data,FUN=sum)
  
  cog_data$cog_depth<-cog_data$dens_depth/cog_data$cpue_kgkm2
  
  cog_data$year_type<-ifelse(cog_data$year %in% c(2016,2018,2019,2022), "warm",'cold')
  cog_data$year_type<-as.factor(cog_data$year_type)                
  
  #only slope years
  cog_data1<-subset(cog_data,year %in% c(2002, 2004, 2008, 2010, 2012, 2016)) 
  cog_data1$year<-as.factor(cog_data1$year)
  cog_data1$year<-droplevels(cog_data1$year)
  
  #cog_data$depth_bin<-factor(cog_data$depth_bin, levels = c("0-50", "50-100", '100-150',"150-200",">200"))
  cog_data$year<-as.factor(cog_data$year)
  
  p1<-
    ggplot()+
    geom_point(data=cog_data1,aes(x=year,y=cog_depth,group=year),position=position_dodge(width=0.5))+
    #scale_color_brewer(name = 'depth bins (m)')+
    theme_bw()+
    theme(axis.title.x = element_blank(),axis.text.x = element_text(color = ifelse(levels(cog_data1$year) %in% c(2016,2018,2019,2022), "red",
                                                                                   ifelse(levels(cog_data1$year) %in% c(1992,2007:2010,2012), "blue",'black')),size=14),plot.title = element_text(size=16),legend.text = element_text(size = 12),  # Increase legend text size
          legend.title = element_text(size = 14), # Increase legend title size
          legend.key.size = unit(1.5, 'lines'),axis.title.y = element_text(size=14),
          legend.position = 'none')+
    ggtitle(label =s)+
    ylab(label = 'depth COG (m)')#+  
  
  p2<-
    ggplot()+
      geom_boxplot(data=cog_data1,aes(x=year_type,y=cog_depth,group=year_type),position=position_dodge(width=0.8,preserve = "single"),outliers = FALSE)+
      #scale_fill_brewer(name = 'depth bins (m)')+
      theme_bw()+
      theme(axis.title.x = element_blank(),axis.text.x = element_text(color = ifelse(levels(cog_data1$year_type) == 'warm', "red",'blue'),size=14),axis.title.y = element_blank(),plot.title = element_text(size=16),legend.text = element_text(size = 12),  # Increase legend text size
            legend.title = element_text(size = 14), # Increase legend title size
            legend.key.size = unit(1.5, 'lines'))+
      ggtitle(label ='')#+
  
  # print(
  #   cowplot::plot_grid(p1,p2,rel_widths = c(1,0.6))
  # )
  

  #####################
  # MAX DEPTH
  #####################
  
  data_nonzero<-subset(data,cpue_kgkm2!=0)
  
  max_data<-
    aggregate(depth_m ~ year ,data_nonzero,FUN=max)
  
  max_data$year_type<-ifelse(max_data$year %in% c(2016,2018,2019,2022), "warm",'cold')
  max_data$year_type<-as.factor(max_data$year_type)                
  
  #only slope years
  max_data1<-subset(max_data,year %in% c(2002, 2004, 2008, 2010, 2012, 2016)) 
  max_data1$year<-as.factor(max_data1$year)
  max_data1$year<-droplevels(max_data1$year)
  
  #cog_data$depth_bin<-factor(cog_data$depth_bin, levels = c("0-50", "50-100", '100-150',"150-200",">200"))
  max_data1$year<-as.factor(max_data1$year)
  
  p1<-
    ggplot()+
    geom_point(data=max_data1,aes(x=year,y=depth_m,group=year),position=position_dodge(width=0.5))+
    #scale_color_brewer(name = 'depth bins (m)')+
    theme_bw()+
    theme(axis.title.x = element_blank(),axis.text.x = element_text(color = ifelse(levels(max_data1$year) %in% c(2016,2018,2019,2022), "red",
                                                                                   ifelse(levels(max_data1$year) %in% c(1992,2007:2010,2012), "blue",'black')),size=14),plot.title = element_text(size=16),legend.text = element_text(size = 12),  # Increase legend text size
          legend.title = element_text(size = 14), # Increase legend title size
          legend.key.size = unit(1.5, 'lines'),axis.title.y = element_text(size=14),
          legend.position = 'none')+
    ggtitle(label =s)+
    ylab(label = 'depth COG (m)')#+  
  
  p2<-
    ggplot()+
    geom_boxplot(data=max_data1,aes(x=year_type,y=depth_m,group=year_type),position=position_dodge(width=0.8,preserve = "single"),outliers = FALSE)+
    #scale_fill_brewer(name = 'depth bins (m)')+
    theme_bw()+
    theme(axis.title.x = element_blank(),axis.text.x = element_text(color = ifelse(levels(max_data1$year_type) == 'warm', "red",'blue'),size=14),axis.title.y = element_blank(),plot.title = element_text(size=16),legend.text = element_text(size = 12),  # Increase legend text size
          legend.title = element_text(size = 14), # Increase legend title size
          legend.key.size = unit(1.5, 'lines'))+
    ggtitle(label ='')#+
  
  # print(
  #   cowplot::plot_grid(p1,p2,rel_widths = c(1,0.6))
  # )
  
 }
#load(file = './data processed/grid_EBS_NBS.RData')


################
# EFFECTIVE AREA OCCUPIED
################

#load grid of NBS and EBS
load('./extrapolation grids/northern_bering_sea_grid.rda')
load('./extrapolation grids/eastern_bering_sea_grid.rda')
grid<-as.data.frame(rbind(data.frame(northern_bering_sea_grid,region='NBS'),data.frame(eastern_bering_sea_grid,region='EBS')))
grid$cell<-1:nrow(grid)

#HIST simulated data
load(file = paste0('./output/species/ms_sim_dens.RData'))  #sim_dens1

for (s in sp_shelf) {
  
  s<-sp_shelfslope[1]
  
  #####################
  # DENSITY DATA
  #####################
  
  sim<-1
  
  ind<-colSums(sim_dens1[,s,,sim]*grid$Area_in_survey_km2)
  ave_dens<-colMeans(sim_dens1[,s,,sim])
  
  area<-ind/ave_dens
  area1<-reshape2::melt(as.matrix(area))
  names(area1)<-c('year','dummy','area')
  
  area1<-area1[which(area1$year %in% c(1992,2007:2010,2012,2016,2018,2019,2022)),]
  
  area1$year<-as.factor(area1$year)
  
  
  ggplot()+
    geom_point(data = area1,aes(x=year,y=area))+
    theme_bw()+
    theme(axis.title.x = element_blank(),axis.text.x = element_text(color = ifelse(levels(area1$year) %in% c(2016,2018,2019,2022), "red",
                                                                                   ifelse(levels(area1$year) %in% c(1992,2007:2010,2012), "blue",'black')),size=14),plot.title = element_text(size=16),legend.text = element_text(size = 12))  # Increase legend text size
          
    
  ####
  #### check depths plots with loess
  ####
  #### 
  
  
}





################
# lower bound depth percentile
################

################
# cog slope by using temp ebs
################

load(file='./output/cog_OM.RData') #cog_df

#average SBT by region and year
temp_region<-
  aggregate(Temp ~ Year + region,grid.ebs_year,FUN=mean)

names(temp_region)<-c('year','region','temp')

temp_region1<-subset(temp_region,region=='EBSshelf')
temp_region1$region<-'EBSslope'

###################
# merge 
###################
cog_df1<-subset(cog_df,region=='EBSslope')
cog_df2<-merge(cog_df1,temp_region,by=c('year','region'))

#plot

cog_df2$coloryear<-ifelse(cog_df2$year %in%  c(2016,2018,2019,2022), "red",
                          ifelse(cog_df2$year %in% c(1992,2007:2010,2012),'blue','grey'))

###################
# plot region temp vs depth cog over species and regions 
###################

# Define a function to scale values
# scale_value <- function(x) {
#   (x - mean(x)) / sd(x)
# }
# 
# 
# # Add a scaled column using group_by and mutate
# cog_df11 <- cog_df1 %>%
#   group_by(sp, region) %>%
#   mutate(scale_bio = scale_value(bio))

#p<-
  ggplot() +
  geom_point(data = cog_df2, aes(x = cog_depth, y = temp, color = as.factor(year)),alpha=0.6) +
  facet_wrap(region ~ sp, scales = 'free', nrow = 2) +
  labs(x = 'bathymetric center of gravity', y = 'mean region SBT') +
  scale_color_manual(values = unique(cog_df2[,c("coloryear","year")])[,'coloryear'], 
                     breaks = unique(cog_df2[,c("coloryear","year")])[,'year'],  # Specify the breaks based on unique year values
                     labels = unique(cog_df2[,c("coloryear","year")])[,'year'],  # Use unique year values as labels
                     name = 'year') +
  #geom_vline(data=depth_region,aes(xintercept=as.numeric(DepthGEBCO)),linetype='dotted')+
  theme_bw()+
  theme(legend.position = 'none')+
  #geom_smooth(data = cog_df1, aes(x = depth, y = temp), method = "lm", se = FALSE, color = 'black',alpha=0.9)
  geom_line(data = cog_df2, aes(x = cog_depth, y = temp),stat="smooth",method = "lm", 
            linetype ="dashed",
            color='black')+
  scale_size() #limits = c(100000,1000000000000),




#########################
# RANGE EDGES BY DEPTH AREAS
#########################

#grid without slope
grid.ebs_year1shelf<-subset(grid.ebs_year1,region %in% c('EBSshelf','NBS'))

#loop over common spp in shelf and slope
for (s in sp_shelfslope) {
  
  #s<-sp_shelfslope[1]#;y<-colnames(dens_slope)[1]
  
  
  #EBS+NBS model fit and density
  load(paste0('./shelf EBS NBS VAST/',s,'/fit.RData'))
  fit_shelf<-fit
  bio_shelf<-fit_shelf$Report$Index_gctl[,,,1]
  
  plot_list<-list()
  
  #loop over years
  #warm and cold years #c(1992,2007:2010,2012,2016,2018,2019,2022)
  wcyrs<-c(1992,2007:2010,2012,2016,2018,2019,2022)
  
  for (y in as.factor(c(1992,2007:2010,2012,2016,2018,2019,2022))) {
    
    cat(paste('#################### ',s," - ",y,'\n'))
    
    
    #grid without slope
    grid.ebs_year1shelf<-subset(grid.ebs_year,region %in% c('EBSshelf','NBS') & Year==y)
    
    all_bio2<-c(bio_shelf[,y])
    all_bio3<-cbind(all_bio2,grid.ebs_year1shelf)
    all_bio3$dens<-all_bio3$all_bio2/all_bio3$Area_in_survey_km2
    all_bio4<-subset(all_bio3,region=='EBSshelf')
    
    p5<-quantile(all_bio4$dens,c(0.10))
    p95<-quantile(all_bio4$dens,c(0.90))
    
    
    all_bio4$quantile<-ifelse(all_bio4$dens < p95 & all_bio4$dens > p5,1,0)
    all_bio4$all<-1
    
    #df to spatialpoint df
    coordinates(all_bio4) <- ~ Lon + Lat
    crs(all_bio4)<-c(crs="+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
    
    #reproject coordinates for plotting purposes
    all_bio5<-spTransform(all_bio4,'+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs')
    all_bio6<-data.frame(all_bio5)
    
    #x and y cells
    xycells<-as.integer(sqrt(dim(all_bio5)[1]))
    
    # create a template raster
    r1 <- raster(ext=extent(all_bio5),ncol=xycells, nrow=xycells) #c(15800,15800) 7000
    
    #create raster
    r2<-rasterize(all_bio5, r1 ,field=c('quantile','Temp','all'))
    crs(r2) <- CRS('+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs')
    r2$quantile[r2$quantile==0]<-NA
    r2$Temp[r2$Temp>2]<-NA
    r2$Temp[r2$Temp<2]<-1
    
    
    #create polygon to get boundaries of each strata
    r3<-as.data.frame(r2,xy=TRUE)
    r3<-r3[complete.cases(r3$quantile),] 
    #as factors
    r3$quantile<-as.factor(r3$quantile)
    r4a<-rasterToPolygons(r2$quantile,dissolve=TRUE,digits = 1)
    r4b<-rasterToPolygons(r2$Temp,dissolve=TRUE,digits = 1)
    r4c<-rasterToPolygons(r2$all,dissolve=TRUE,digits = 1)
    
    #plot
    p<-
      ggplot()+
      geom_polygon(data=r4b,aes(x=long,y=lat,group=group,color='cold pool'), fill='blue',alpha=0.7)+
      geom_raster(data=r3,aes(x=x,y=y,fill=quantile),alpha=0.4)+
      scale_fill_manual(values='grey30',labels='species range')+
      geom_polygon(data=r4a,aes(x=long,y=lat,group=group),fill='transparent', colour="grey50")+
      scale_color_manual(values=c('cold pool'='transparent'))+
      #geom_polygon(data=ak_sppoly,aes(x=long,y=lat,group=group),color='transparent',linewidth=0.2,fill = 'bisque3')+
      geom_polygon(data=EBSshelf_sh,aes(x=long,y=lat,group=group),fill=NA,col='black')+
      theme_bw()+
      labs(title = y)+
      annotate("text", x = -500000, y = 1500000, label = y, size = 6) +
      theme(axis.title = element_blank(),
            legend.title = element_blank(),
            plot.margin = unit(c(0, 0, 0, 0),units = 'in'), 
            axis.text = element_blank(),
            legend.position = 'none',
            panel.grid = element_blank(),
            plot.title = element_blank(),
            axis.ticks = element_blank())
    
    plot_list[[y]]<-p
    
    
    
    
    
    
    
    
    
    
    # #loop over regions
    # for (r in c(na.omit(unique(grid.ebs_year1$region_depth)),'200_400','0_200')) {
    #   
    #   if (r=='0_100') {
    #     c<-which(grid.ebs_year1$region_depth=="0_100")
    #   } else if (r=='100_200') {
    #     c<-which(grid.ebs_year1$region_depth=='100_200')
    #   } else if (r == '200_300') {
    #     c <- which(grid.ebs_year1$region_depth=='200_300')
    #   } else if (r=="300_400") {
    #     c<-which(grid.ebs_year1$region_depth=='300_400')
    #   } else if (r=='200_400') {
    #     c<-which(grid.ebs_year1$region_depth %in% c('200_300','300_400'))
    #   } else if (r=='0_200') {
    #     c<-which(grid.ebs_year1$region_depth %in% c('0_100','100_200'))
    #   }
    #   
    #   #subset grid and dens by grid
    #   igrid<-grid.ebs_year1[c,]
    #   all_bio3<-all_bio2[c]
    #   
    #   #calculate percentiles
    #   #data.frame(t(quantile(all_bio3, probs = c(0.05, 0.10, 0.90, 0.95)))) # Calculate percentiles
    #   
    #   #bio to dens
    #   all_dens3<-all_bio2[c]/igrid$Area_in_survey_km2      
    #   
    #   #append
    #   ind_df<-rbind(ind_df,
    #                 data.frame('sp'=s,
    #                            'year'=y,
    #                            'abundance'=sum(all_bio3),
    #                            'dens'=mean(all_dens3),
    #                            'region'=r,
    #                            'percentiles'=data.frame(t(quantile(all_bio3, probs = c(0.05, 0.10, 0.90, 0.95)))))) # Calculate percentiles)))  
    
    
    
  }
  
  
  legend1 <- cowplot::get_legend( 
    plot_list[[1]] + 
      theme(legend.position = "right",legend.spacing.y = unit(-0.1,'in')) #,legend.justification = c(0.72,0.5)
  ) 
  
  #plot(legend1)
  
  
  plot_list[['2023']]<-legend1
  
  agg_png(paste0('./figures slope/',s,'_rangeEBS_filtered.png'), width = 10, height = 10, units = "in", res = 300)
  print(
    plot_grid(plotlist = plot_list)
  )
  dev.off()
  
}
  





for (s in sp_shelfslope) {
  
  s<-sp_shelfslope[2]
  
  #####################
  # DENSITY DATA
  #####################
  
  #shelf
  data<-readRDS(paste0('./data processed/species/',s,'/data_geostat_temp.rds'))
  
  #filter by region and years
  data_shelf1<-data[which(data$survey_name=='Eastern Bering Sea Crab/Groundfish Bottom Trawl Survey'),]
  data_shelf2<-data_shelf1[which(data_shelf1$year %in% c(1992,2007:2010,2012,2016,2018,2019,2022)),]
  data_shelf3<-data_shelf2[,c("lat_start","lon_start","year", "cpue_kgkm2" , "depth_m")]
  
  #slope
  data_slope1<-data[which(data$survey_name=='Eastern Bering Sea Slope Bottom Trawl Survey'),]
  data_slope2<-data_slope1[which(data_slope1$year %in% c(1992,2007:2010,2012,2016,2018,2019,2022)),]
  data_slope3<-data_slope2[,c("lat_start","lon_start","year", "cpue_kgkm2" , "depth_m")]
  #EBS_slope: 2002, 2004, 2008, 2010, 2012 and 2016 
  
  data<-rbind(data_shelf3,data_slope3)
  
  data$year_type<-ifelse(data$year %in% c(2016,2018,2019,2022), "warm",'cold')
  data$year_type<-as.factor(data$year_type)                
  data$year<-as.factor(data$year)
  
  data1<-subset(data,cpue_kgkm2!=0)
  
  # Load necessary libraries
  library(dplyr)
  library(ggplot2)
  library(tidyr)
  
  mean_depths <- data1 %>%
    group_by(year_type) %>%
    summarise(MeanDepth = mean(depth_m, na.rm = TRUE)) %>%
    pivot_wider(names_from = year_type, values_from = MeanDepth)
  
  # Calculate median depths for each species in warm and cold years
  median_depths <- data1 %>%
    group_by(year_type) %>%
    summarise(MedianDepth = median(depth_m, na.rm = TRUE)) %>%
    pivot_wider(names_from = year_type, values_from = MedianDepth)
  
  
  print(paste0(s))
  
  print(
  median_depths
  )
  
  
  # Kernel density estimation for depth distributions
  ggplot(data, aes(x = depth_m, fill = year_type)) +
    geom_density(alpha = 0.5) +
    #facet_wrap(~ Species) +
    labs(x = 'Depth (m)', y = 'Density', title = 'Depth Density Distribution in Warm vs. Cold Years') +
    scale_fill_manual(values = c('warm' = 'red', 'cold' = 'blue'), name = 'Year Category') +
    theme_minimal()
  

}



all_data<-data.frame(matrix(ncol = 7,nrow = 0))
names(all_data)<-c('lat_start' ,'lon_start' ,'year'  ,'cpue_kgkm2', 'depth_m', 'year_type' ,'species')

for (s in sp_shelfslope) {
  
  #s<-sp_shelfslope[2]
  
  #####################
  # DENSITY DATA
  #####################
  
  #shelf
  data<-readRDS(paste0('./data processed/species/',s,'/data_geostat_temp.rds'))
  
  #filter by region and years
  data_shelf1<-data[which(data$survey_name=='Eastern Bering Sea Crab/Groundfish Bottom Trawl Survey'),]
  data_shelf2<-data_shelf1[which(data_shelf1$year %in% c(1992,2007:2010,2012,2016,2018,2019,2022)),]
  data_shelf3<-data_shelf2[,c("lat_start","lon_start","year", "cpue_kgkm2" , "depth_m")]
  
  #slope
  data_slope1<-data[which(data$survey_name=='Eastern Bering Sea Slope Bottom Trawl Survey'),]
  data_slope2<-data_slope1[which(data_slope1$year %in% c(1992,2007:2010,2012,2016,2018,2019,2022)),]
  data_slope3<-data_slope2[,c("lat_start","lon_start","year", "cpue_kgkm2" , "depth_m")]
  #EBS_slope: 2002, 2004, 2008, 2010, 2012 and 2016 
  
  data<-rbind(data_shelf3,data_slope3)
  
  data$year_type<-ifelse(data$year %in% c(2016,2018,2019,2022), "warm",'cold')
  data$year_type<-as.factor(data$year_type)                
  data$year<-as.factor(data$year)
  data$species<-s
  
  all_data<-rbind(all_data,data)
  
}

summary(all_data)

all_data1<-subset(all_data,cpue_kgkm2!=0)

# Create a new column for weighted depths
all_data1 <- all_data %>%
  mutate(WeightedDepth = depth_m * cpue_kgkm2)

all_data1<-subset(all_data1,depth_m<800)

# Kernel density estimation for depth distributions
ggplot(all_data1, aes(x = depth_m, fill = year_type)) +
  geom_density(alpha = 0.5) +
  facet_wrap(~ species,scales='free_y') +
  labs(x = 'Depth (m)', y = 'Density', title = 'Depth Density Distribution in Warm vs. Cold Years') +
  scale_fill_manual(values = c('warm' = 'red', 'cold' = 'blue'), name = 'Year Category') +
  theme_minimal()



