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
pack_cran<-c('googledrive','lubridate')

#install pacman to use p_load function - call library and if not installed, then install
if (!('pacman' %in% installed.packages())) {
  install.packages("pacman")}

#load/install packages
pacman::p_load(pack_cran,character.only = TRUE)

#setwd
out_dir<-'C:/Users/Daniel.Vilas/Work/Adapting Monitoring to a Changing Seascape/'  
out_dir<-'/Users/daniel/Work/Adapting Monitoring to a Changing Seascape/'
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
           'Paralithodes camtschaticus')

#get files from google drive and set up
files<-googledrive::drive_find()
32 #for dvilasg@uw.edu

#get id shared folder from google drive
id.bering.folder<-files[which(files$name=='Bering redesign RWP project'),'id']

#list of files and folder
files.1<-googledrive::drive_ls(id.bering.folder$id)
id.data<-files.1[which(files.1$name=='data raw'),'id']
files.2<-googledrive::drive_ls(id.data$id)

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

#filter by species
catch1<-subset(catch,scientific_name %in% spp)
length(unique(catch1$scientific_name))==length(spp)

#####################################
# MERGE CATCH and HAUL DATA 
#####################################
#if there are 15 selected spp and 16693 hauls in slope EBS, shelf EBS and NBS
#then 15*16693=250395 rows for the dataframe

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

#save data_geostat file
saveRDS(all1, paste0('./data processed/species/slope_shelf_EBS_NBS_data_geostat.rds'))

#loop over species to create data_geostat df
for (sp in spp) {
  
  #sp<-spp[3]
  
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



