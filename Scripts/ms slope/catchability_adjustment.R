####################################################################
####################################################################
##
##    Extract sea bottom temperature (SBT) from netcdf of Bering 10K ROMS
##    Netcdf downloaded from https://data.pmel.noaa.gov/aclim/thredds/
##    Incorporate Temp (SBT) to Bering Sea grid
##    Incorporate Temp (SBT) to data_geostat file from Bering Sea 
##    (shelf Eastern Bering Sea, slope Eastern Bering Sea, northern Bering Sea)
##    Save data_geostat_temp file to fit VAST
##    Daniel Vilas (danielvilasgonzalez@gmail.com/dvilasg@uw.edu/daniel.vilas@noaa.gov)
##    Lewis Barnett, Stan Kotwicki, Zack Oyafuso, Megsie Siple, Leah Zacher, Lukas Defilippo, Andre Punt
##
####################################################################
####################################################################

#clear all objects
rm(list = ls(all.names = TRUE)) 
#free up memrory and report the memory usage
gc() 

#libraries from cran to call or install/load
pack_cran<-c('ncdf4','raster','FNN','lubridate','ggpubr')

#install pacman to use p_load function - call library and if not installed, then install
if (!('pacman' %in% installed.packages())) {
  install.packages("pacman")}

#load/install packages
pacman::p_load(pack_cran,character.only = TRUE)

#setwd - depends on computer using
#out_dir<-'C:/Users/Daniel.Vilas/Work/Adapting Monitoring to a Changing Seascape/' #NOAA laptop  
out_dir<-'/Users/daniel/Work/Adapting Monitoring to a Changing Seascape/' #mac
#out_dir<-'/Users/daniel/Work/VM' #VM
setwd(out_dir)

#range years of data
# sta_y<-1982
# end_y<-2022

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
       #'Lepidopsetta sp.',
       'Chionoecetes bairdi',
       'Sebastes alutus',
       'Sebastes melanostictus',
       'Atheresthes evermanni',
       'Sebastes borealis',
       'Sebastolobus alascanus',
       'Glyptocephalus zachirus',
       'Bathyraja aleutica')

#remove Anoploma and Reinhardtius because habitat preference reasons
#spp<-setdiff(spp, c('Anoplopoma fimbria','Reinhardtius hippoglossoides'))

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
        'Saffron cod',
        'Sablefish',
        'Snow crab',
        'Blue king crab',
        'Red king crab',
        'Tanner crab',
        'Pacific ocean perch',
        'Rougheye and blackspotted rockfish',
        'Kamchatka flounder',
        'Shortraker rockfish',
        'Shortspine thornyhead',
        'Rex sole',
        'Aleutian skate')

#read length raw data
data_length<-readRDS('./data raw/ak_bts_ebs_nbs_slope.rds') #data_length
head(data_length)
head(data_length$specimen)
dim(data_length$specimen)

#get cruisejoin for the ebs
head(data_length$cruise)
cruisejoin_ebs<-subset(data_length$cruise,SURVEY=='EBS')[,'CRUISEJOIN']

#get hauls for the ebs
hauls_ebs<-subset(data_length$haul,CRUISEJOIN %in% cruisejoin_ebs)

#convert time
time_axis <- as.POSIXct(hauls_ebs$START_TIME, origin = "1900-01-01", tz = "GMT") 
hauls_ebs$YEAR <- format(time_axis, "%Y")

#check hauls number per year in ebs
aggregate(HAULJOIN ~ YEAR,hauls_ebs,FUN=length)

#code species
spp_code<-unique(data_length$species[,c('SPECIES_CODE',"REPORT_NAME_SCIENTIFIC")])
names(spp_code)<-c('species_code',"scientific_name")
spp_code1<-spp_code[which(spp_code$scientific_name %in% 
                            spp),]

#data input 
coef_wl<-expand.grid('spp'=spp_code1$scientific_name,
            'sex'=c('1','2','all'),
            'year'=as.character(c(2002,2004,2008,2010,2012,2016)))

#to add values parms
coef_wl$log_a<-NA
coef_wl$b<-NA

for (r in 1:nrow(coef_wl)) {
  
  #r<-1
  
  #species
  sp<-coef_wl[r,'spp']
  sp_code<-spp_code1[which(spp_code1$scientific_name==sp),'species_code']
  
  #sex
  sex<-coef_wl[r,'sex']
  sex<-ifelse(sex=='1',1,
              ifelse(sex=='2',2,c(1,2)))
  
  #year
  y<-coef_wl[r,'year']
  
  #hauls per year
  hauls<-subset(hauls_ebs,YEAR %in% y)[,'HAULJOIN']
  
  #data weight-length
  data<-subset(data_length$specimen,
               SPECIES_CODE %in% sp_code  & 
                 HAULJOIN %in% hauls  &
                 SEX %in% sex)
  
  #remove NAs and zeros
  data<-data[complete.cases(data),]
  data<-subset(data,LENGTH != 0 & WEIGHT!=0)
  
  #jump if no more than 3 obs
  if (nrow(data)<3) {
    next
  }
  
  #fit lm log space
  m <- stats::lm(log(WEIGHT) ~ log(LENGTH), data = data)
  pars <- as.list(coef(m))
  pars <- stats::setNames(pars, c("log_a", "b"))
  
  #store parms
  coef_wl[r,'log_a']<-pars$log_a
  coef_wl[r,'b']<- pars$b 

}

#remove species with no coeffs
sp_data<-unique(coef_wl[complete.cases(coef_wl),'spp'])
sp_rem<-setdiff(spp_code1$scientific_name,sp_data)
  
#get cruisejoin for the slope
head(data_length$cruise)
cruisejoin_bss<-subset(data_length$cruise,SURVEY=='BSS')[,'CRUISEJOIN']

#get hauls for the ebs
hauls_bss<-subset(data_length$haul,CRUISEJOIN %in% cruisejoin_bss)

#convert time
time_axis <- as.POSIXct(hauls_bss$START_TIME, origin = "1900-01-01", tz = "GMT") 
hauls_bss$YEAR <- format(time_axis, "%Y")

#length data slope
length_bss<-data_length$size[which(data_length$size$HAULJOIN %in% unique(hauls_bss$HAULJOIN)),]
length_bss<-merge(length_bss,hauls_bss[,c('CRUISEJOIN','HAULJOIN','YEAR')],by=c('CRUISEJOIN','HAULJOIN'))
length_bss<-merge(length_bss,spp_code1,by.x='SPECIES_CODE',by.y='species_code')

#only for species with coeffs
length_bss<-length_bss[which(length_bss$scientific_name %in% sp_data),]

#add weight
length_bss$WEIGHT<-NA

for (r in 1:nrow(length_bss)) {
  
  #r<-1
  cat(paste('#',r,'-',nrow(length_bss),'\n'))
  
  #length
  l<-length_bss[r,'LENGTH']
  
  #year
  y<-length_bss[r,'YEAR']
  
  #sex
  sex<-as.character(length_bss[r,'SEX'])
  sex<-ifelse(sex=='3','all',sex)
  
  #spp
  sp<-length_bss[r,'scientific_name']
  
  #filter
  coef_wl1<-coef_wl[which(coef_wl$spp==sp & coef_wl$sex==sex & coef_wl$year==y),]
  
  #if no data, get the average
  if (is.na(coef_wl1$log_a)) {
    coef_wl1<-coef_wl[which(coef_wl$spp==sp & coef_wl$sex==sex ),]
    coef_wl1<-data.frame('spp'=sp,'sex'=sex,'year'=NA,'log_a'=mean(coef_wl1$log_a,na.rm=TRUE),'b'=mean(coef_wl1$b,na.rm=TRUE))
  }
  
  #convert length to weight  
  length_bss[r,'WEIGHT'] <- exp(coef_wl1$log_a + coef_wl1$b * log(l))
}

#weight of all individuals at length bin
length_bss$WEIGHT_FREQ<-length_bss$WEIGHT*length_bss$FREQUENCY

#apply catchability ratio at length, sex and year to obtain the weight as it would be in the EBS
2

#weigth by species for each haul
wl<-aggregate(WEIGHT_FREQ ~ scientific_name + YEAR + HAULJOIN,length_bss,FUN=sum)

#check data_catch slope
catch_bss<-data_length$catch[which(data_length$catch$HAULJOIN %in% unique(hauls_bss$HAULJOIN)),]
catch_bss<-merge(catch_bss,hauls_bss[,c('HAULJOIN','YEAR')],by=c('HAULJOIN'))
catch_bss<-merge(catch_bss,spp_code1,by.x='SPECIES_CODE',by.y='species_code')

#weigth by species for each haul
wc<-aggregate(WEIGHT ~ scientific_name + YEAR + HAULJOIN,catch_bss,FUN=sum)

#merge both to check how similar
merge(wc,wl,by=c('scientific_name','YEAR','HAULJOIN'))

