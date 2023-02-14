####################################################################
####################################################################
##
##    get raw data from the EBS and NBS 
##    (shelf Eastern Bering Sea, slope Eastern Bering Sea, northern Bering Sea)
##    add sp names
##    include depth from GEBCO
##    create data_geostat file to fit VAST
##    Daniel Vilas (danielvilasgonzalez@gmail.com/dvilasg@uw.edu)
##
####################################################################
####################################################################

#clear all objects
rm(list = ls(all.names = TRUE)) 
#free up memrory and report the memory usage
gc() 

#libraries from cran to call or install/load
pack_cran<-c('googledrive','raster')

#install pacman to use p_load function - call library and if not installed, then install
if (!('pacman' %in% installed.packages())) {
  install.packages("pacman")}

#load/install packages
pacman::p_load(pack_cran,character.only = TRUE)

#setwd
out_dir<-'E:/UW/Adapting Monitoring to a Changing Seascape/'
setwd(out_dir)

#range years of data
sta_y<-2002
end_y<-2022

#get files from google drive and set up
files<-googledrive::drive_find()
1 #for dvilasg@uw.edu

#get id shared folder from google drive
id.bering.folder<-files[which(files$name=='Bering redesign RWP project'),'id']

#list of files and folder
files.1<-googledrive::drive_ls(id.bering.folder$id)
id.data<-files.1[which(files.1$name=='data raw'),'id']
files.2<-googledrive::drive_ls(id.data$id)

#####################################
# SP CODE
#####################################

#create directory
dir.create('./Data/Additional/',showWarnings = FALSE)

#get species file
files.bering<-files.2$name
file<-files.bering[grep('species',files.bering)]
file.id<-files.2[which(files.2$name==file),'id']

#download file into temp folder
googledrive::drive_download(file=file.id$id,
                            path = paste0('./Data/Additional/',file),
                            overwrite = TRUE)

#read csv file
df.sp<-read.csv(paste0('./Data/Additional/',file))
df.sp1<-df.sp[,c('SPECIES_CODE','SPECIES_NAME')]

#sp list
sp.list<-c('Limanda aspera','Gadus chalcogrammus','Gadus macrocephalus','Atheresthes stomias','Reinhardtius hippoglossoides',
           'Lepidopsetta polyxystra','Hippoglossoides elassodon','Pleuronectes quadrituberculatus','Hippoglossoides robustus')
           #'Boreogadus saida','Eleginus gracilis','Anoplopoma fimbria',
           #'Chionoecetes opilio','Paralithodes platypus','Paralithodes camtschaticus'

#####################################
# CPUE DATA
#####################################

#create directory
dir.create('./Data/Surveys/',showWarnings = FALSE)

#get survey folder id
id.bering.folder<-files[which(files$name=='Surveys'),'id']

#list of files and folder
files.1<-googledrive::drive_ls(id.bering.folder$id)

#get cpue file
file<-files.1[grep('|slope|shelf|nbs',files.1),]
#file.id<-files.2[which(files.2$name %in% file),]

#df to store results
dfsurveys<-data.frame(matrix(nrow = 0,ncol = 6))
colnames(dfsurveys)<-c("SPECIES_NAME","YEAR","LATITUDE","LONGITUDE","CPUE_KGHA",'SURVEY')

 for (i in 1:nrow(file)) {
   
   #i=2
   
   googledrive::drive_download(file=file$id[i],
                               path = paste0('./Data/Surveys/',file$name[i]),
                               overwrite = TRUE)
   
   
   df<-read.csv(paste0('./Data/Surveys/',file$name[i]))
   
   if (file$name[i]=="ebs_slope_cpue.csv") {
     
     #add sp information
     df<-merge(x = df,
               y = df.sp1,
               by = 'SPECIES_CODE',
               all.x = TRUE)
     
     #rename column
     names(df)[10]<-'CPUE_KGHA'
     
     #correct slope data (slope is in ha, while others and in km2)
     df$CPUE_KGHA<-df$CPUE_KGHA/100}
   
   #filter by sp
   df1<-subset(df, SPECIES_NAME %in% sp.list) 
   
   #select columns
   df2<-df1[,c("SPECIES_NAME","YEAR","LATITUDE","LONGITUDE","CPUE_KGHA")]
   
   #add survey column
   
   if (file$name[i]=='nbs_cpue.csv') {
     df2$SURVEY<-'nbs'
   } else if (file$name[i]=='ebs_slope_cpue.csv') {
     df2$SURVEY<-'slope'
   } else if (file$name[i]=='ebs_shelf_cpue.csv') {
     df2$SURVEY<-'shelf'}
   
  #rbind data 
  dfsurveys<-rbind(dfsurveys,df2)
   
 }

#####################################
# ADD DEPTH
#####################################

#create directory
dir.create('./Data/Bathymetry/',showWarnings = FALSE)

#get raster depth from gebco (https://download.gebco.net/)
#get id shared folder from google drive
id.bering.folder<-files[which(files$name=='Bathymetry'),'id']

#list of files and folder
files.1<-googledrive::drive_ls(id.bering.folder$id)
id.data<-files.1[which(files.1$name=='gebco_2022_n70.0_s50.0_w-180.0_e-155.0.asc'),]

googledrive::drive_download(file=id.data$id,
                            path = paste0('./Data/Bathymetry/',id.data$name),
                            overwrite = TRUE)

#read raster
r<-raster('./Data/Bathymetry/gebco_2022_n70.0_s50.0_w-180.0_e-155.0.asc')

#extract depth values for each station
rr<-extract(r, SpatialPoints(cbind(dfsurveys$LONGITUDE,dfsurveys$LATITUDE)))
dfsurveys$DEPTH<--rr

#rename cols
colnames(dfsurveys)<-c('Species','Year',"Lat","Lon",'CPUE_kg','Survey','Depth')

#scale grid bathymetry values to standard normal, using the mean and sd
dfsurveys$LogDepth <- log(dfsurveys$Depth)
dfsurveys$ScaleLogDepth <- scale(dfsurveys$LogDepth)

#####################################
# CREATE DATA_GEOSTAT FILE
#####################################

#create folder
dir.create('./slope shelf EBS NBS VAST/',showWarnings = FALSE)

#save data_geostat file
saveRDS(dfsurveys, paste0('./slope shelf EBS NBS VAST/slope_shelf_ebs_nbs_data_geostat.rds'))

#loop over species to create data_geostat df
for (sp in sp.list) {
  
  #sp<-sp.list[3]
  
  #print species to check progress
  cat(paste("    -----", sp, "-----\n"))
  
  #create folder to store results
  dir.create(paste0('./slope shelf EBS NBS VAST/',sp),
             showWarnings = FALSE)
  
  #filter by sp
  df1<-subset(dfsurveys, Species == sp)
  df1<-subset(df1, Year %in% sta_y:end_y)
  
  #remove rows with NAs in env data
  df1<-df1[complete.cases(df1[c('Depth')]),]

  #save data_geostat file
  saveRDS(df1, paste0('./slope shelf EBS NBS VAST/',sp,'/data_geostat.rds'))
  
}

