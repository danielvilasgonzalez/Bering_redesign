####################################################################
####################################################################
##
##    extract sea bottom temperature (SBT) from netcdf of Bering 10K ROMS
##    incorporate SBT to Bering Sea grid
##    incorporate SBT to data_geostat file from Bering Sea 
##    (shelf Eastern Bering Sea, slope Eastern Bering Sea, northern Bering Sea)
##    save data_geostat_temp file to fit VAST
##    Daniel Vilas (danielvilasgonzalez@gmail.com/dvilasg@uw.edu)
##
####################################################################
####################################################################

#clear all objects
rm(list = ls(all.names = TRUE)) 
#free up memrory and report the memory usage
gc() 

#libraries from cran to call or install/load
pack_cran<-c('ncdf4','raster','FNN')

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

#####################################
# BERING SEA GRIDS
#####################################
 
#create directory
dir.create('./Data/Extrapolation Grids/',showWarnings = FALSE)
 
#get files from google drive and set up
files<-googledrive::drive_find()
1 #for dvilasg@uw.edu

#get id shared folder from google drive
id.bering.folder<-files[which(files$name=='Extrapolation Grids'),'id']

#grid names
grfiles<-c('SlopeThorsonGrid.csv','EBSThorsonGrid.csv','NBSThorsonGrid.csv')

#list of files and folder
files.1<-googledrive::drive_ls(id.bering.folder$id)
id.data<-files.1[which(files.1$name %in% grfiles),]

  #loop over grid files
  for (i in 1:nrow(id.data)) {
    
    #i=2
   googledrive::drive_download(file=id.data$id[i],
                               path = paste0('./Data/Extrapolation Grids/',id.data$name[i]),
                               overwrite = TRUE)
    
   gr<-read.csv(paste0('./Data/Extrapolation Grids/',id.data$name[i]))
   colnames(gr)[10]<-'Area_km2'
   
   #shapefile name
   if (id.data$name[i]=='SlopeThorsonGrid.csv') {
     grname<-'EBSslope_gr'
   } else if (id.data$name[i]=='EBSThorsonGrid.csv') {
     grname<-'EBSshelf_gr'     
   } else if (id.data$name[i]=='NBSThorsonGrid.csv') {
     grname<-'NBS_gr' 
   }
  
   #name strata
   gr$STRATA<-gsub('_gr','',grname)
   
   #assign shapefiles
   assign(grname,gr)
   
  }
 
#join grids
grid.ebs<-rbind(EBSshelf_gr,EBSslope_gr,NBS_gr)
