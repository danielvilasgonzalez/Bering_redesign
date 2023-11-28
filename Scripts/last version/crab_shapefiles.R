library(sf)

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
summary(gdb_table1)
plot(gdb_table1)

gdb_layers$name

for (a in gdb_layers$name) {
  
  
  
  
}


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

#set working directory
out_dir<-'C:/Users/Daniel.Vilas/Work//Adapting Monitoring to a Changing Seascape/'
out_dir<-'/Users/daniel/Work/Adapting Monitoring to a Changing Seascape/'
setwd(out_dir)

# Define plot extent (through trial end error) units km
panel_extent <- data.frame(x = c(-1716559.21, -77636.05), #x = c(-1326559.21, -87636.05),
                           y = c(483099.5, 2194909.7)) #y = c(533099.5, 1894909.7))

#Alaska layer
ebs_layers <- akgfmaps::get_base_layers(select.region = "ebs", set.crs = "+proj=longlat +datum=WGS84 +no_defs")
ak_sppoly<-as(ebs_layers$akland, 'Spatial')



gdb_layers$name

gdb_table1


ggplot()+
  geom_sf(data=ebs_layers$survey.strata,fill = NA)+ #gdb_table
  # geom_sf(data=gdb_table1,fill = 'green',alpha=0.1)+
  # geom_sf_text(data = gdb_table1, label=gdb_layers$name[1],color='green')+
  # # geom_sf(data=gdb_table2,fill = 'red',alpha=0.1)+
  # # geom_sf_text(data = gdb_table2, label=gdb_layers$name[2],color='red')+
  # geom_sf(data=gdb_table3,fill = 'yellow',alpha=0.1)+
  # geom_sf_text(data = gdb_table3, label=gdb_layers$name[3],color='yellow')+
  # # geom_sf(data=gdb_table4,fill = 'blue',alpha=0.1)+
  # # geom_sf_text(data = gdb_table4, label=gdb_layers$name[4],color='blue')+
  # geom_sf(data=gdb_table5,fill = 'orange',alpha=0.1)+
  # geom_sf_text(data = gdb_table5, label=gdb_layers$name[5],color='orange')+
  geom_sf(data=gdb_table6,fill = 'black',alpha=0.1)+
  geom_sf_text(data = gdb_table6, label=gdb_layers$name[6],color='black')+
  #geom_polygon(data=ak_sppoly,aes(x=long,y=lat,group=group),fill = 'grey60')+
  scale_x_continuous(expand = c(0,0),position = 'bottom',
                     breaks = c(-175,-170,-165,-160,-155),sec.axis = dup_axis())+
  #geom_polygon(data=NBS_sh,aes(x=long,y=lat,group=group),fill=NA,col='black')+
  #geom_polygon(data=EBSshelf_sh,aes(x=long,y=lat,group=group),fill=NA,col='black')+
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













