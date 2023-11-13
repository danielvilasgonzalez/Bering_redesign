####################################################################
####################################################################
##
##    get raw bottom trawl (BT) data from the EBS and NBS 
##    (shelf Eastern Bering Sea, slope Eastern Bering Sea, northern Bering Sea)
##    add sp names
##    include depth from GEBCO
##    create data_geostat file to fit OM VAST
##    Daniel Vilas (danielvilasgonzalez@gmail.com/dvilasg@uw.edu)
##
####################################################################
####################################################################

#clear all objects
rm(list = ls(all.names = TRUE)) 
#free up memrory and report the memory usage
gc() 

#libraries from cran to call or install/load
pack_cran<-c('googledrive','lubridate','ggplot2','fishualize')

#install pacman to use p_load function - call library and if not installed, then install
if (!('pacman' %in% installed.packages())) {
  install.packages("pacman")}

#load/install packages
pacman::p_load(pack_cran,character.only = TRUE)

#install akgfmaps to extract shapefile of Alaska
if (!('akgfmaps' %in% installed.packages())) {
  devtools::install_github('afsc-gap-products/akgfmaps')};library(akgfmaps)

#setwd - depends on computer using
#out_dir<-'C:/Users/Daniel.Vilas/Work/Adapting Monitoring to a Changing Seascape/' #NOAA laptop  
#out_dir<-'/Users/daniel/Work/Adapting Monitoring to a Changing Seascape/' #mac
out_dir<-'/Users/daniel/Work/VM' #VM
setwd(out_dir)

#range years of data
sta_y<-1982
end_y<-2022

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
           'Paralithodes camtschaticus',
           'Lepidopsetta sp.',
           'Chionoecetes bairdi',
           'Sebastes alutus',
           'Sebastes melanostictus')

#get files from google drive and set up
files<-googledrive::drive_find()
3 #2 #for dvilasg@uw.edu

#get id shared folder from google drive
id.bering.folder<-files[which(files$name=='Bering redesign RWP project'),'id']

#list of files and folder
files.1<-googledrive::drive_ls(id.bering.folder$id)
id.data<-files.1[which(files.1$name=='data raw'),'id']
files.2<-googledrive::drive_ls(id.data$id)

#create directory
dir.create('./extrapolation grids/',showWarnings = FALSE)

#get id shared folder from google drive
id.bering.folder<-files[which(files$name=='Bering redesign DV'),'id']

#list of files and folder
files.1<-googledrive::drive_ls(id.bering.folder$id)
id.data<-files.1[which(files.1$name=='extrapolation grids'),'id']
files.2<-googledrive::drive_ls(id.data$id)

#download file
#eastern
googledrive::drive_download(file=files.2$id[3],
                            path = paste0('./extrapolation grids/',files.2$name[3]),
                            overwrite = TRUE)
#northern
googledrive::drive_download(file=files.2$id[4],
                            path = paste0('./extrapolation grids/',files.2$name[4]),
                            overwrite = TRUE)
#slope
googledrive::drive_download(file=files.2$id[5],
                            path = paste0('./extrapolation grids/',files.2$name[5]),
                            overwrite = TRUE)

#####################################
# HAUL DATA
#####################################

#create directory
dir.create('./data raw/',showWarnings = FALSE)

#get haul (stations) data
file<-files.2[grep('haul',files.2$name),]
#file.id<-files.2[which(files.2$name %in% file),]

#download file
googledrive::drive_download(file=file$id,
                            path = paste0('./data raw/',file$name),
                            overwrite = TRUE)

#read csv file
haul<-readRDS(paste0('./data raw/',file$name))
dim(haul);length(unique(haul$hauljoin))

#####################################
# CATCH DATA
#####################################

#get cpue file
file<-files.2[grep('catch',files.2$name),]
#file.id<-files.2[which(files.2$name %in% file),]

#download file
googledrive::drive_download(file=file$id,
                            path = paste0('./data raw/',file$name),
                            overwrite = TRUE)

#read csv file
catch<-readRDS(paste0('./data raw/',file$name))

#check for names - based on important spp for the slope (all that we have + Blackspotted/rougheye and POP)
unique(catch$common_name)[grep('perch',unique(catch$common_name))] #Pacific ocean perch
unique(catch$common_name)[grep('blackspotted',unique(catch$common_name))] #"rougheye and blackspotted rockfish unid." "blackspotted rockfish"
unique(catch[which(catch$common_name=='Pacific ocean perch'),'scientific_name']) #Sebastes alutus
unique(catch[which(catch$common_name=='rougheye and blackspotted rockfish unid.'),'scientific_name']) #NA - to arrange
unique(catch[which(catch$common_name=='blackspotted rockfish'),'scientific_name']) #Sebastes melanostictus

#add scientific_name to 'rougheye and blackspotted rockfish unid.'
catch$scientific_name[catch$common_name == 'rougheye and blackspotted rockfish unid.'] <- 'Sebastes melanostictus'

#most northern rock sole was missidentified before 1996
unique(catch$common_name)[grepl('rock sole',unique(catch$common_name))]
subset(catch, common_name=='rock sole unid.')

#filter by species
catch1<-subset(catch,scientific_name %in% spp)

#sum blackspotted rockfish and blackspotted rockfish unid
catch2<-catch1[which(catch1$scientific_name=='Sebastes melanostictus'),]
catch21<-aggregate(catch2[, c('cpue_kgha','cpue_kgkm2','cpue_noha','cpue_nokm2','count','weight_kg')], 
                   by = list('hauljoin'=catch2$hauljoin), FUN = sum)
catch3<-cbind('hauljoin'=catch21$hauljoin,
              'species_code'=unique(catch[which(catch$common_name=='blackspotted rockfish'),'species_code']),
              catch21[,-1],
              'taxon_confidence'='Unassessed',
              'scientific_name'='Sebastes melanostictus',
              'common_name'='rougheye and blackspotted rockfish',
              'worms'=unique(catch[which(catch$common_name=='blackspotted rockfish'),'worms']),
              'itis'=NA)
catch1<-subset(catch1,scientific_name != "Sebastes melanostictus")
catch1<-rbind(catch1,catch3)

length(unique(catch1$scientific_name))==length(spp)

#####################################
# MERGE CATCH and HAUL DATA 
#####################################
#if there are 19 selected spp and 16693 hauls in slope EBS, shelf EBS and NBS
#then 19*16693=317167 rows for the dataframe

#create the empty df 
haul1<-do.call("rbind", replicate(length(spp), haul, simplify = FALSE))
dim(haul1)

#replicate spp for each station
spp1<-rep(spp,each=nrow(haul))

#join dataframe
all<-data.frame(haul1,'scientific_name'=spp1)
head(all);dim(all)

#merge haul and catch
all1<-merge(all,catch1,all.x=T)
dim(all1)

#cpue_kgha,cpue_kgkm2,cpue_noha,cpue_nokm2,count,weight_kg columns need to replace NA by 0s
all1[c('cpue_kgha','cpue_kgkm2','cpue_noha','cpue_nokm2','count','weight_kg')][is.na(all1[c('cpue_kgha','cpue_kgkm2','cpue_noha','cpue_nokm2','count','weight_kg')])] <- 0
head(all1)
summary(all1)

#####################################
# CREATE DATA_GEOSTAT FILE
#####################################

#create folder
dir.create('./data processed/',showWarnings = FALSE)
dir.create('./data processed/species/',showWarnings = FALSE)

#add year and month
all1$month<-month(as.POSIXlt(all1$date, format="%d/%m/%Y"))
all1$year<-year(as.POSIXlt(all1$date, format="%d/%m/%Y"))

#check Lepidopsetta sp. 
mm<-subset(all1,scientific_name=='Lepidopsetta sp.')
mm$year<-lubridate::year(mm$date)
tapply(mm$count,mm$year,summary)

#remove Lepidopsetta sp. >=1996
all1<-all1[-which(all1$scientific_name=='Lepidopsetta sp.' & all1$year>=1996),]

#remove Lepidopsetta sp. <=1995
all1<-all1[-which(all1$scientific_name=='Lepidopsetta polyxystra' & all1$year<=1995),]

#replace spp rock sole unid
all1$scientific_name[all1$scientific_name == 'Lepidopsetta sp.'] <- 'Lepidopsetta polyxystra'

#remove from spp vector
spp<-spp[spp!='Lepidopsetta sp.']

#save data_geostat file
saveRDS(all1, paste0('./data processed/species/slope_shelf_EBS_NBS_data_geostat.rds'))

#loop over species to create data_geostat df
for (sp in spp) {

  #sp<-spp[16]
  
  #print species to check progress
  cat(paste("    -----", sp, "-----\n"))
  
  #create folder to store results
  dir.create(paste0('./data processed/species/',sp),
             showWarnings = FALSE)
  
  #filter by sp
  all2<-subset(all1, scientific_name == sp)
  all2<-subset(all2, year %in% sta_y:end_y)
  cat(paste("    ----- ", nrow(all2) , "samples -----\n"))
  
  #xx<-all2[which(is.na(all2$bottom_temp_c)),]
  #summary(xx)
  #save data_geostat file
  saveRDS(all2, paste0('./data processed/species/',sp,'/data_geostat.rds'))
  
}

#####################################
# PLOT 
#####################################

haul$year<-year(as.POSIXlt(haul$date, format="%d/%m/%Y"))

#mean SBT by year
yearagg.df <- aggregate(data = haul, bottom_temp_c ~ year, mean)
yearagg.df$TempAnomaly<-NA

#calculate SBT anomaly
for (i in 2:nrow(yearagg.df)) {
  yearagg.df[i,'TempAnomaly']<-yearagg.df[i,'bottom_temp_c']-yearagg.df[i-1,'bottom_temp_c']
}

#plot 1
print(
  ggplot() +
    geom_line(data=yearagg.df, aes(x=year, y=bottom_temp_c),linetype='dashed')+
    geom_point(data=yearagg.df, aes(x=year, y=bottom_temp_c,color=bottom_temp_c),size=2)+
    geom_bar(data=yearagg.df, aes(x=year, y=TempAnomaly,fill=TempAnomaly),color='black',stat="identity",position = position_dodge(0.9))+
    scale_colour_gradient2(low = 'darkblue',high='darkred',midpoint = 2.5)+
    scale_fill_gradient2(low = 'darkblue',high='darkred',midpoint = 0)+
    #xlab(label = 1982:2022)+
    scale_x_continuous(breaks=c(1982:2022),expand = c(0,0.1))+
    theme_bw()+
    labs(y='°C',color='SBT',fill='SBT anomaly',x='')+
    guides(fill=guide_legend(order = 2),color=guide_legend(order = 1))+
    theme(panel.grid.minor.x = element_blank(),legend.spacing  = unit(1,'cm'),
          panel.grid.minor.y = element_line(linetype=2,color='grey90'),axis.text.x = element_text(angle=90,vjust=0.5))
)

#plot 2
print(
  ggplot() +
    geom_density(data=haul, aes(x=bottom_temp_c, group=year))+
    geom_vline(data=yearagg.df,aes(xintercept=bottom_temp_c,group=year, color=bottom_temp_c),
               linetype="dashed", size=1)+
    scale_x_continuous(limits = c(-2,8),breaks =c(-1,1,3,5,7))+
    scale_colour_gradient2(low = 'darkblue',high='darkred',midpoint = 2.5)+
    #geom_text(data=yearagg.df,aes(label=paste0('mean BotTemp = ',round(bottom_temp_c,digits = 2)),x = 10,y=0.43))+
    facet_wrap(~year,nrow = 3)+
    labs(x='°C',y='',color='SBT')+
    theme_bw())
#}



#####################################
# PLOT SPECIES ESTIMATES PER YEAR AND REGION
#####################################

head(all1)

#add common name
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
       'Paralithodes camtschaticus',
       'Chionoecetes bairdi',
       'Sebastes alutus',
       'Sebastes melanostictus')

#remove Anoploma and Reinhardtius because habitat preference reasons
#spp<-setdiff(spp, c('Anoplopoma fimbria','Reinhardtius hippoglossoides'))

#remove two species because of habitat preferece reasons
all1<-all1[all1$scientific_name %in% spp,]

#common names
spp1<-c('Yellowfin sole',
        'Alaska pollock',
        'Pacific cod',
        'Arrowtooth flounder',
        'Greenland turbot',
        'Northern rock sole',
        'Flathead sole',
        'Alaska plaice',
        'Bering flounder',
        'Arctic cod',
        'Saffon cod',
        'Sablefish',
        'Snow crab',
        'Blue king crab',
        'Red king crab',
        'Tanner crab',
        'Pacific ocean perch',
        'Rougheye and blackspotted rockfish')

#df sp scientific and common
df_spp<-data.frame('spp'=spp,
                   'common'=spp1)

#merge both df
all2<-merge(all1,df_spp,by.x='scientific_name',by.y = 'spp',all.x = 'TRUE')

#get sci + common name
all2$sp<-paste0(all2$scientific_name,'\n(',all2$common,')')
all2$scientific_name2<-gsub(' ','_',all2$scientific_name)

#merge lepidosettas
all1$scientific_name[all1$scientific_name == 'Lepidopsetta sp.'] <- 'Lepidopsetta polyxystra'

#plot CPUE
p<-
ggplot()+
  geom_boxplot(data=all2,aes(x=year,y=cpue_kgha,group=interaction(year,survey_name),color=survey_name),alpha=0.7)+
  facet_wrap(~sp,scales = 'free_y',nrow=5)+
  #add_fishape(data=all2,aes(option = scientific_name2))+
  scale_color_manual(values=c("Northern Bering Sea Crab/Groundfish Survey - Eastern Bering Sea Shelf Survey Extension"="#4682B4",
                              "Eastern Bering Sea Crab/Groundfish Bottom Trawl Survey"="#B4464B",
                              "Eastern Bering Sea Slope Bottom Trawl Survey"="#B4AF46"),
                     labels = c('EBS shelf','EBS slope','NBS'),name='survey')+
  scale_x_continuous(breaks = c(1985,1990,1995,2000,2005,2010,2015,2020),
                     minor_breaks = setdiff(1982:2022,c(1982,1985,1990,1995,2000,2005,2010,2015,2020,2022)))+
  scale_y_continuous(limits = c(0,NA),labels = scales::comma)+
  theme_bw()+
  labs(y='CPUE',x='')+
  theme(panel.grid.minor = element_line(linetype=2,color='grey90'),strip.background = element_rect(fill='white'),
        legend.position=c(.80,.08),legend.key.size = unit(20, 'points'),legend.text = element_text(size=10),
        legend.title = element_text(size=14),strip.text = element_text(size=12))+ #axis.text.x = element_text(angle=90,vjust=0.5),
  expand_limits(y = 0)

#save plot
ragg::agg_png(paste0('./figures/CPUE_survey_year.png'), width = 13, height = 10, units = "in", res = 300)
p
dev.off()


#plot CPUE in log+1
p<-
  ggplot()+
  geom_boxplot(data=all2,aes(x=year,y=log(cpue_kgha+1),group=interaction(year,survey_name),color=survey_name),alpha=0.7,position = position_dodge2(preserve = "single"))+
  facet_wrap(~sp,scales = 'free_y',nrow=5)+
  #add_fishape(data=all2,aes(option = scientific_name2))+
  scale_color_manual(values=c("Northern Bering Sea Crab/Groundfish Survey - Eastern Bering Sea Shelf Survey Extension"="#4682B4",
                              "Eastern Bering Sea Crab/Groundfish Bottom Trawl Survey"="#B4464B",
                              "Eastern Bering Sea Slope Bottom Trawl Survey"="#B4AF46"),
                     labels = c('EBS shelf','EBS slope','NBS'),name='survey')+
  scale_x_continuous(breaks = c(1985,1990,1995,2000,2005,2010,2015,2020),
                     minor_breaks = setdiff(1982:2022,c(1982,1985,1990,1995,2000,2005,2010,2015,2020,2022)))+
  scale_y_continuous(limits = c(0,NA),labels = scales::comma)+
  theme_bw()+
  labs(y='log(CPUE+1)',x='')+
  theme(panel.grid.minor = element_line(linetype=2,color='grey90'),strip.background = element_rect(fill='white'),
        legend.position=c(.80,.08),legend.key.size = unit(20, 'points'),legend.text = element_text(size=10),
        legend.title = element_text(size=14),strip.text = element_text(size=12))+ #axis.text.x = element_text(angle=90,vjust=0.5),
  expand_limits(y = 0)

#save plot
ragg::agg_png(paste0('./figures/CPUE_survey_year.png'), width = 13, height = 10, units = "in", res = 300)
p
dev.off()

#################################################
# CREATE DATA SAMPLING SCENARIO BASELINE
#################################################

#load grid of NBS and EBS
load('./extrapolation grids/northern_bering_sea_grid.rda')
load('./extrapolation grids/eastern_bering_sea_grid.rda')
grid<-as.data.frame(rbind(data.frame(northern_bering_sea_grid,region='NBS'),data.frame(eastern_bering_sea_grid,region='EBS')))
grid$cell<-1:nrow(grid)
#df to spatialpoint df
coordinates(grid) <- ~ Lon + Lat
crs(grid)<-c(crs='+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs')
#reproject coordinates for plotting purposes
D2_1<-grid
D2_2<-data.frame(D2_1)
#x and y cells
xycells<-as.integer(sqrt(dim(D2_1)[1]))
# create a template raster
r1 <- raster(ext=extent(D2_1),ncol=xycells, nrow=xycells) #c(15800,15800) 7000
#create raster
r2<-rasterize(D2_1, r1 ,field='cell')
#plot(r2)

#EBS and NBS layer
ebs_layers <- akgfmaps::get_base_layers(select.region = "ebs", set.crs = "EPSG:3338")

#baseline strata areas
strata_areas<-as.data.frame(ebs_layers$survey.strata)
#dataframe stratum and area
strata_areas <- data.frame('Stratum'=strata_areas$Stratum,'Area_in_survey_km2'=strata_areas$Precise_Ar/1000)

#strata polygon
strata_pol<-as(ebs_layers$survey.strata, 'Spatial')
proj4string(strata_pol) <- CRS('+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs') 
strata_pol<-spTransform(strata_pol,CRSobj = CRS("+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs "))

#locations
x<-ebs_layers$survey.grid
st<-x$STATIONID
xx<-st_geometry(x)
xxx<-st_centroid(xx)
#plot(xxx);class(xxx)
coords <- st_coordinates(xxx)
lat <- coords[, 2]
lon <- coords[, 1]
# plot(lon,lat)
# text(lon, lat, st, pos = 3)
baseline<-data.frame('Lat'=lat,'Lon'=lon,'stationid'=st)
#corner stations
corner<- c('GF','HG','IH','QP','JI','ON','PO')
st.corner<-paste(corner,collapse = '|')
baseline$corner<-ifelse(grepl(st.corner,baseline$stationid),TRUE,FALSE)
#locations of stations
locations <- as.data.frame(baseline)
st<-baseline
coordinates(st)<- ~ Lon + Lat
proj4string(st) <- CRS('+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs') 
st<-spTransform(st,CRSobj = CRS('+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs'))
#cell
locations$cell<-extract(r2,st)
st1<-as.data.frame(st)[,c("Lon","Lat")]
names(st1)<-c('x','y')
xy<-st1
sampled = apply(X = xy, MARGIN = 1, FUN = function(xy) r2@data@values[which.min(replace(distanceFromPoints(r2,xy), is.na(r2), NA))])
locations$cell<-sampled
locations$Stratum<-over(st,strata_pol)[,'Stratum']

#number of samples per strata for random sampling
y<-aggregate(locations$cell,by=list(locations$Stratum),length)
yc<-aggregate(subset(locations,corner!=TRUE)[,'cell'],by=list(subset(locations,corner!=TRUE)[,'Stratum']),length)
n_samples<-data.frame('stratum'=yc$Group.1,'scnbase'=y$x,'scnbase_bis'=yc$x)

# grid1<-grid
# coordinates(grid1)<- ~ Lon + Lat
# crs(grid1)<-'+proj=aea +lat_0=50 +lon_0=-154 +lat_1=55 +lat_2=65 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs'
# 
# xx<- as(x, 'Spatial')
# #grid over polygon to get samples grid which current baseline strata
# cell_strata<-data.frame(as.data.frame(grid1,'stratum'=over(grid1,xx)[,'Stratum']))

#list baseline strata
baseline_strata<-list(strata_areas=strata_areas,locations=locations,n_samples=n_samples,cell_strata=as.data.frame(grid))

#create directory
dir.create('./output/',showWarnings = FALSE)
#save data
save(baseline_strata,file='./output/baseline_strata.RData')
