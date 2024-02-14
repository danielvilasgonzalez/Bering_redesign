library(sf)
library(sp)

# Set the path to the geodatabase file
gdb_path <- "./shapefiles/CrabStrataShapefiles_GAPgrid.gdb/"

# List the layers/tables in the geodatabase file
gdb_layers <- st_layers(gdb_path)

# Select the specific table you want to read
selected_layer <- gdb_layers$name[1]

# Read the selected table from the geodatabase
#"BBRKC_strata"        "Pribilof_BKC_strata" "Pribilof_RKC_strata" "StMatt_BKC_strata"   "Norton_RKC_Strata"  
#[6] "EBS_CO_CB_strata"   
gdb_table1 <- st_read(dsn = gdb_path, layer = gdb_layers$name[1])
gdb_table2 <- st_read(dsn = gdb_path, layer = gdb_layers$name[2])
gdb_table3 <- st_read(dsn = gdb_path, layer = gdb_layers$name[3])
gdb_table4 <- st_read(dsn = gdb_path, layer = gdb_layers$name[4])
gdb_table5 <- st_read(dsn = gdb_path, layer = gdb_layers$name[5])
gdb_table6 <- st_read(dsn = gdb_path, layer = gdb_layers$name[6])

# View the structure and summary of the data
str(gdb_table1)
crs(gdb_table1)

summary(gdb_table1)
plot(gdb_table1$Shape)
gdb_table1$Shape


#load grid of NBS and EBS
load('./extrapolation grids/northern_bering_sea_grid.rda')
load('./extrapolation grids/eastern_bering_sea_grid.rda')
grid<-as.data.frame(rbind(data.frame(northern_bering_sea_grid,region='NBS'),data.frame(eastern_bering_sea_grid,region='EBS')))
grid$cell<-1:nrow(grid)
grid2<-grid

area<-list()

for (i in 1:length(gdb_layers$name)) {
  
  #i<-3
  
  #load grid of NBS and EBS
  load('./extrapolation grids/northern_bering_sea_grid.rda')
  load('./extrapolation grids/eastern_bering_sea_grid.rda')
  grid<-as.data.frame(rbind(data.frame(northern_bering_sea_grid,region='NBS'),data.frame(eastern_bering_sea_grid,region='EBS')))
  grid$cell<-1:nrow(grid)
  #grid2<-grid
  coordinates(grid)<-~Lon + Lat
  proj4string(grid) <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0") 
  pts<-spTransform(grid,CRSobj = crs(gdb_table1))
  #pts<-as.data.frame(pts)
  

  #st to sp
  gdb_tablea<-as(get(paste0('gdb_table',i)),'Spatial')
  #plot(pts)
  plot(gdb_tablea)
  
  if (length(gdb_tablea$Shape_Area)!=1) {
    iarea<-sum(gdb_tablea$Shape_Area/1000000)
  } else{
    iarea<-gdb_tablea$Shape_Area/1000000
  }
  
  area[[i]]<-iarea
    
  
  plot(gdb_tablea[2])
  
  #points over polygon
  xx<-over(pts,gdb_tablea)
  
  # Check which points fall within the sf object using st_within
  points_within_polygon <- as.data.frame(pts)[!is.na(xx), 'cell']
  
  #name<-
  grid2$newcolumn<-FALSE
  grid2$newcolumn[points_within_polygon]<-TRUE
  names(grid2)[ncol(grid2)]<-gdb_layers$name[i]
}

area[[1]]<-NULL
area[[4]]<-area[[5]]

grid2

ggplot()+
  geom_point(data=grid2,aes(x=Lon,y=Lat,col=EBS_CO_CB_strata))


samp<-'scn1'
#load survey allocations by sampling design
load(file = paste0('./output/survey_allocations_',samp,'.RData')) #scn_allocations
dimnames(scn_allocations)[[3]]<-c('sys','rand','sb')
dimnames(scn_allocations)


PBL_KC_cells<-grid2[which(grid2$Pribilof_BKC_strata==TRUE),'cell']
STM_BKC_cells<-grid2[which(grid2$StMatt_BKC_strata==TRUE),'cell']
EBS_C_cells<-grid2[which(grid2$EBS_CO_CB_strata==TRUE),'cell']

sur<-1
apr<-'sys'
sim<-1
y<-1

#spp-stock to add
crabs<-c('PBL_BKC','PBL_RKC','STM_BKC','SNW_CRB','TNR_CRB')

#load data
load(file = paste0('./output/species/ms_sim_dens.RData'))  #sim_dens1
#simulated densities of survey and year
sim_dens2<-sim_dens1[,,1,sim]

scn_allocations[scn_allocations[,'sur',apr]==sur & scn_allocations[,'cell',apr] %in% n,,apr]




#spp-stock to add
crabs<-c('PBL_BKC','PBL_RKC','STM_BKC','SNW_CRB','TNR_CRB')
names(area)<-crabs
crabs_spp<-c('Paralithodes platypus','Paralithodes camtschaticus','Paralithodes platypus','Chionoecetes opilio','Chionoecetes bairdi')


for (c in crabs) {
  
  c<-crabs[1]
  
  if (grepl('PBL',c)) {
    cells<-PBL_KC_cells
  } else if (grepl('CRB',c)) {
    cells<-EBS_C_cells
  } else {
    cells<-STM_BKC_cells
  }
  
  crab_sp<-crabs_spp[match(c,crabs)]
  
  #get densities based on station allocations
  sim_survey_crabs<-data.frame(cbind(strata=scn_allocations[scn_allocations[,'sur',apr]==n & scn_allocations[,'cell',apr] %in% cells,c('strata'),apr],
                               dens=sim_dens2[scn_allocations[scn_allocations[,'sur',apr]==n & scn_allocations[,'cell',apr] %in% cells,c('cell'),apr],crab_sp]),check.names = FALSE)
  names(sim_survey_crabs)[ncol(sim_survey_crabs)]<-crab_sp
  
  
  sim_survey1<-reshape2::melt(sim_survey_crabs,id.vars=c('strata'))
  sim_survey1$strata<-'crab'
  
  #mean, sum and var by strata and year (variable)
  sim_survey2<-aggregate(x=sim_survey1$value,
                         by=list(strata=sim_survey1$strata,sp=sim_survey1$variable),
                         FUN = function(x) c('mean' = mean(x,na.rm=T), 'sum' = sum(x),'var' = var(x,na.rm=T) ))
  zzz<-data.frame('strata'=sim_survey2$strata,'sp'=sim_survey2$sp,'mean'=sim_survey2$x[,c('mean')],'var'=sim_survey2$x[,c('var')]) #/length(yy$value)
  
  #add index strata for sum to compute index (mean strata density * area of strata) kg!
  zzz$index_strata<-zzz$mean*area[[c]]
  
  #add strata var 
  zzz$strs_var<-zzz$var*(area[[c]]^2)/nrow(sim_survey_crabs) #sum(survey_detail$Nh) 
  
  #get CV across years
  zzz$cv<- sqrt(zzz$strs_var) / zzz$index_strata
  
  #mean CV 
  mean(zzzz1$cv,na.rm=TRUE)
  
  #get outputs
  STRS_mean <- zzzz1$index_strata
  STRS_var <- zzzz1$strs_var
  CV <- sqrt(STRS_var) / STRS_mean
  
}



#get densities based on station allocations
sim_survey<-data.frame(cbind(strata=scn_allocations[scn_allocations[,'sur',apr]==n,c('strata'),apr],
                             dens=sim_dens2[scn_allocations[scn_allocations[,'sur',apr]==n,c('cell'),apr],]),check.names = FALSE)

sim_survey1<-reshape2::melt(sim_survey,id.vars=c('strata'))
#mean, sum and var by strata and year (variable)
sim_survey2<-aggregate(x=sim_survey1$value,
                       by=list(strata=sim_survey1$strata,sp=sim_survey1$variable),
                       FUN = function(x) c('mean' = mean(x,na.rm=T), 'sum' = sum(x),'var' = var(x,na.rm=T) ))




print(selected_rows)

####################################################################
####################################################################
##
##    Create a Bering Sea maps
##    Sampling approach maps
##    Daniel Vilas (danielvilasgonzalez@gmail.com/dvilasg@uw.edu)
##
####################################################################
####################################################################

#####################################
# Settings
#####################################


#libraries from cran to call or install/load
pack_cran<-c('ggspatial','raster','rasterVis','rgeos','scales','rnaturalearth','grid','ggplot2')

#install pacman to use p_load function - call library and if not installed, then install
if (!('pacman' %in% installed.packages())) {
  install.packages("pacman")}

#install akgfmaps to extract shapefile of Alaska
if (!('akgfmaps' %in% installed.packages())) {
  devtools::install_github('afsc-gap-products/akgfmaps')};library(akgfmaps)

#load/install packages
pacman::p_load(pack_cran,character.only = TRUE)
out_dir<-'/Users/daniel/Work/Adapting Monitoring to a Changing Seascape/'# Define plot extent (through trial end error) units k
sur<-1
apr<-'sys'

sur<-1
apr<-'sys'
m
panel_extent <- data.frame(x = c(-1716559.21, -77636.05), #x = c(-1326559.21, -87636.05),
                           y = c(483099.5, 2194909.7)) #y = c(533099.5, 1894909.7))

#Alaska layer
#Alaska layer
ebs_layers <- akgfmaps::get_base_layers(select.region = "ebs", set.crs = "EPSG:3338")
ak_sppoly<-as(ebs_layers$akland, 'Spatial')


p<-
ggplot()+
  #geom_sf(data=ebs_layers$survey.strata,fill = NA)+ #gdb_table
  geom_polygon(data=ak_sppoly,aes(x=long,y=lat,group=group),color='black',linewidth=0.2,fill = 'grey80')+
  geom_sf(data=gdb_table1,fill = 'green',alpha=0.1)+
  geom_sf_text(data = gdb_table1, label='Bristol\nBay',color='green',nudge_y = 20000,lineheight = 0.9)+
  geom_sf(data=st_union(gdb_table2),fill = 'red',alpha=0.1)+
  geom_sf_text(data = st_union(gdb_table2), label='Pribilof\nIslands',color='red',nudge_y = 70000,lineheight = 0.9)+
  #geom_sf(data=st_union(gdb_table3),fill = 'yellow',alpha=0.1)+
  #geom_sf_text(data = st_union(gdb_table3), label='Pribilof\nIslands',color='yellow',nudge_x = 60000,nudge_y = 50000,lineheight = 0.9)+
  geom_sf(data=st_union(gdb_table4),fill = 'blue',alpha=0.1)+
  geom_sf_text(data = st_union(gdb_table4), label='St. Matthew\nIsland',lineheight = 0.9,color='blue',nudge_x = 50000,nudge_y = -30000)+
  #geom_sf(data=gdb_table5,fill = 'orange',alpha=0.1)+
  #geom_sf_text(data = gdb_table5, label=gdb_layers$name[5],color='orange')+
  geom_sf(data=st_union(gdb_table6),fill = 'black',alpha=0.1)+
  geom_sf_text(data = st_union(gdb_table6), label='EBS',color='black',nudge_y = 90000)+
  geom_polygon(data=ak_sppoly,aes(x=long,y=lat,group=group),fill = 'grey60')+
  scale_x_continuous(expand = c(0,0),position = 'bottom',
                     breaks = c(-175,-170,-165,-160,-155),sec.axis = dup_axis())+
  geom_polygon(data=NBS_sh,aes(x=long,y=lat,group=group),fill=NA,col='black')+
  geom_polygon(data=EBSshelf_sh,aes(x=long,y=lat,group=group),fill=NA,col='black')+
  coord_sf(crs = '+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs',
           #xlim = c(panel_extent$x[1]+200000,panel_extent$x[2]),
           xlim = c(panel_extent$x[1]+200000,panel_extent$x[2]+100000),
           ylim = c(panel_extent$y[1]-100000,panel_extent$y[2]-200000),
           label_axes = "-NE-")+
  theme(aspect.ratio = 1,panel.grid.major = element_line(color = rgb(0, 0, 0,20, maxColorValue = 285), linetype = 'dashed', linewidth =  0.5),
        panel.background = element_rect(fill = NA),panel.ontop = TRUE,text = element_text(size=10),
        legend.background =  element_rect(fill = "transparent", colour = "transparent"),legend.key.height= unit(25, 'points'),
        legend.key.width= unit(25, 'points'),axis.title = element_blank(),legend.position = c(0.12,0.47), #c(0.12,0.47)
        panel.border = element_rect(fill = NA, colour = 'black'),legend.key = element_rect(color="black"),
        axis.text = element_text(color='black'),legend.spacing.y = unit(10, 'points'),
        axis.text.y.right = element_text(hjust= 0.1 ,margin = margin(0,7,0,-25, unit = 'points'),color='black'),
        axis.text.x = element_text(vjust = 6, margin = margin(-7,0,7,0, unit = 'points'),color='black'),
        axis.ticks.length = unit(-5,"points"))+
  scale_y_continuous(expand = c(0,0),position = 'right',sec.axis = dup_axis())#+


#save plot
ragg::agg_png(paste0('./figures/crabs_management units.png'), width = 7, height = 7, units = "in", res = 300)
p
dev.off()

########################
# stratification crabs from LDF
#######################

crab_spp<-c('bairdi','bkc','opilio','rkc')

for (sp in crab_spp) {
  
strata<-read.csv(paste0('/Users/daniel/Documents/GitHub/Corner-Stations-II/strata_',sp,'_newts_baseline.csv'))
strata$STRATUM<-as.factor(strata$STRATUM)
# strata20<-data.frame('X',NA,NA,'2020',mean(strata$LATITUDE),mean(strata$LONGITUDE),NA,0)
# names(strata20)<-names(strata)
# strata<-rbind(strata,strata20)
# strata$SURVEY_YEAR<-as.factor(strata$SURVEY_YEAR)              
# levels(strata$SURVEY_YEAR)<-c(levels(strata$SURVEY_YEAR),c('2020'))
ebs_layers$survey.strata<-spTransform(ebs_layers$survey.strata,CRSobj = "+proj=longlat +datum=WGS84 +no_defs")

tail(strata)
print(
ggplot()+
  geom_sf(data=ebs_layers$survey.strata,fill = NA)+
  geom_point(data=subset(strata,SURVEY_YEAR %in% c(1992,2021)),aes(x=LONGITUDE,y=LATITUDE,color=STRATUM),size=0.5)+
  facet_wrap(~SURVEY_YEAR)+
  xlim(c(-180,-155))+
  ylim(c(50,67))+
  # coord_sf(crs = '+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs',
  #          #xlim = c(panel_extent$x[1]+200000,panel_extent$x[2]),
  #          xlim = c(panel_extent$x[1]+200000,panel_extent$x[2]+100000),
  #          ylim = c(panel_extent$y[1]-100000,panel_extent$y[2]-200000),
  #          label_axes = "-NE-")+
  ggtitle(sp)
)
}

yr<-unique(strata$SURVEY_YEAR)
sort(as.numeric(names(which(table(strata$STRATUM[strata$SURVEY_YEAR==yr]) > 1))))



for(y in 1:n_year){
  #need to ensure that at least two stations were sampled for a strata to be included in a given year for sd calculations
  strata_yr_list[[i]][[y]] <- sort(as.numeric(names(which(table(strata$STRATUM[strata$SURVEY_YEAR==yr[y]]) > 1))))
  station_strata_yr_list[[i]][[y]] <- replicate(length(strata_yr_list[[i]][[y]]), vector())
  for(g in 1:length(strata_yr_list[[i]][[y]])){
    strata_area_list[[i]][[y]][g] <- unique(strata$TOTAL_AREA[strata$SURVEY_YEAR==yr[y] & strata$STRATUM==strata_yr_list[[i]][[y]][g]])*3.4299
    station_strata_yr_list[[i]][[y]][[g]] <- strata$STATION_ID[strata$SURVEY_YEAR==yr[y] & strata$STRATUM==strata_yr_list[[i]][[y]][g]]
    #print(length(station_strata_yr_list[[i]][[y]][[g]]))
  }
}













