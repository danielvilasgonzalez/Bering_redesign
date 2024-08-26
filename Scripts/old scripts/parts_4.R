########################
########################
# length check separation
########################
########################

conv_tab<-read.csv('./tables/slope_ebsnbs_convspp.csv')
sp_shelfslope<-conv_tab[which(conv_tab$slope=='There is no evidence that the model is not converged' & conv_tab$EBS_NBS=='There is no evidence that the model is not converged'),'spp']
sp_shelf<-conv_tab[which(conv_tab$EBS_NBS=='There is no evidence that the model is not converged'),'spp']
sp_slope<-conv_tab[which(conv_tab$slope=='There is no evidence that the model is not converged' ),'spp']
spp1<-sp_shelf[c(2,3,5,13)]

data_length<-readRDS('./data raw/ak_bts_ebs_nbs_slope.RDS')
head(data_length)
str(data_length)
data_length$size
data_length$catch



data_code<-readRDS('./data raw/afsc_catch_raw_2023_2_21.rds')
spp_code<-unique(data_code[,c('species_code',"scientific_name")])

spp_code1<-spp_code[which(spp_code$scientific_name %in% spp1),]

#data_length$catch


#temp
#average SBT by region and year
temp_region<-
  aggregate(Temp ~ Year + region,grid.ebs_year,FUN=mean)
names(temp_region)<-c('year','region','temp')




#df to store results
metrics_df<-data.frame(matrix(NA,nrow = 0,ncol=5))
names(metrics_df)<-c('cog_lat','cog_lon','cog_depth','years','species')


for (sp in spp1[1:3]) {

  #sp<-spp1[1]
  
  spp_code2<-spp_code1[which(spp_code1$scientific_name==sp),'species_code']
  data_length$specimen
  
  data_size<-subset(data_length$size,SPECIES_CODE==spp_code2)
  data_catch<-subset(data_length$catch,SPECIES_CODE==spp_code2)
  data<-merge(data_size,data_length$haul,by=c('HAULJOIN','CRUISEJOIN'))
  
  #convert time
  time_axis <- as.POSIXct(data$START_TIME, origin = "1900-01-01", tz = "GMT") 
  data$year <- format(time_axis, "%Y")
  data$juv<-ifelse(data$LENGTH>=300,FALSE,TRUE)
  #data1<-subset(data,LENGTH>=300)
  
  sum_len<-aggregate(FREQUENCY ~ juv + year +HAULJOIN,data,FUN=sum)
    
  library(dplyr)
  
  df_percent <- sum_len %>%
    group_by(year,HAULJOIN) %>%
    mutate(
      PERCENTAGE = FREQUENCY / sum(FREQUENCY) * 100
    ) %>%
    ungroup()  # remove grouping
    
    df_percent1<-subset(df_percent,juv==FALSE)

    #sum_len$PERCENTAGE <- sum_len$FREQUENCY / sum(sum_len$FREQUENCY) * 100
    ggplot()+
      geom_bar(data=sum_len,aes(x=LENGTH,y=FREQUENCY),stat='identity')+
      #labs(title=y)+
      facet_wrap(~year)
    

    
    ggplot()+
      geom_bar(data=df_percent,aes(x=LENGTH,y=PERCENTAGE),stat='identity',color='black',fill='black')+
      labs(title=sp)+
      theme_minimal()+
      scale_x_continuous(breaks = c(0,50,100,150,200,250,300,350,400,450,500,550,600,650,700,750,800,850,900,950))+
      theme(strip.text.y = element_text(angle = 0, hjust = 0),axis.text.y = element_blank(),strip.background = element_blank(),axis.title = element_blank(),plot.title=element_text(hjust=0.5))+
      facet_grid(year~.)
    
  
  data_geostat<-readRDS(file = paste0('./data processed/species/',sp,'/data_geostat_envs.rds'))
  sort(unique(data_geostat$hauljoin))
  sort(unique(data$HAULJOIN))
  
  data_geostat1<-merge(df_percent1,data_geostat,by.x=c('HAULJOIN','year'),by.y=c('hauljoin','year'))
  data_geostat1$cpue_kgkm2_adult<-data_geostat1$cpue_kgkm2*data_geostat1$PERCENTAGE/100
  #filter by >300
  
  for (y in 1982:2019) {
      
  data_geostat2<-subset(data_geostat1,survey_name=='Eastern Bering Sea Crab/Groundfish Bottom Trawl Survey')
  data_geostat3<-subset(data_geostat2,year==y)
  dens<-data_geostat3$cpue_kgkm2_adult
  #calculate COG from 
  
  metrics_df<-rbind(metrics_df,
  data.frame(
  'cog_lat'=sum(data_geostat3$cpue_kgkm2_adult*data_geostat3$lat_start,na.rm = TRUE)/sum(data_geostat3$cpue_kgkm2_adult,na.rm = TRUE),
  'cog_lon'=sum(data_geostat3$cpue_kgkm2_adult*data_geostat3$lon_start,na.rm = TRUE)/sum(data_geostat3$cpue_kgkm2_adult,na.rm = TRUE),
  'cog_depth'=sum(data_geostat3$cpue_kgkm2_adult*data_geostat3$depth_m,na.rm = TRUE)/sum(data_geostat3$cpue_kgkm2_adult,na.rm = TRUE),
  'years'=y,
  'species'=sp))
  
  }
  
}
  
  
  # library('corrplot') #package corrplot
  # M<-cor(data_geostat3[,c("cpue_kgkm2_adult","NCaS","NCaO","depth_m","Temp")]) #plot matrix
  # corrplot(M, method = "circle") #plot matrix
  

#temp region
temp_region1<-subset(temp_region,region=='EBSshelf')
metrics_df1<-merge(metrics_df,temp_region1,by.x=c('years'),by.y='year')

#drop units
metrics_df1<-drop_units(metrics_df1)

#subset
metrics_df2<-subset(metrics_df1,species==spp1[3])

# Define a custom color scale function
custom_colors <- colorRampPalette(c("blue", "white", "darkred"))

ggplot()+
  geom_shadowtext(data=metrics_df2,aes(x=cog_depth,y=cog_lat,color=temp,label=years),fontface='bold',bg.r = 0.05)+
  #geom_text(data=cog_prey3,aes(x=cog_depth,y=NCaO,fill=temp,label=year),fontface='bold')+
  #annotate("shadowtext",x=350,y=max(df$cumulative_biomass)/2,bg.color = "black",color = get_temp_color(itemp, temp_min, temp_max),label=y,
  #         bg.r = 0.05,size=12)+
  theme_bw()+
  xlab('x')+
  labs(title = paste0(unique(metrics_df2$species),' (adult)'),x='bathymetrical COG (m)',y='latitudinal COG')+
  theme(legend.position = 'none',aspect.ratio = 1)+
  scale_color_gradientn(colors = custom_colors(100),name = 'SBT (°C)')







################################################################################################
#########################             RMARKDOWN            ############################
################################################################################################


#Get observed data temp, depth and biomass for 

#####################
# OBSERVED DATA
#####################

spp<-list.files('./data processed/species/')
spp1<-spp[c(1:24)]

obs_df<-data.frame(matrix(NA,nrow=0,ncol=7))
colnames(obs_df)<-c('lat_start' ,'lon_start', 'year','cpue_kgkm2', 'depth_m' ,'survey_name','species')

for (s in spp1) {
  
  #s<-spp1[1]

  data<-readRDS(paste0('./data processed/species/',s,'/data_geostat_temp.rds'))


  
#data1<-data[which(data$survey_name %in% c('Eastern Bering Sea Crab/Groundfish Bottom Trawl Survey')),]
#data2<-data1[which(data1$year %in% c(1999,2008,2009,2010,2012,2016,2018,2022)),]
data3<-data[,c("lat_start","lon_start","year", "cpue_kgkm2" , "depth_m",'survey_name')]

# Select rows where any column is NA
rows_with_na <- data3[apply(data3, 1, function(x) any(is.na(x))), ]
print(rows_with_na)

#data3$year_type<-ifelse(data3$year %in% c(2016,2018,2022), "warm",'cold')
#data$year_type<-as.factor(data$year_type)                
#data$year<-as.factor(data$year)

#data4<-subset(data3,cpue_kgkm2!=0)
data3$species<-s
data3<-data3[complete.cases(data3),]

obs_df<-rbind(obs_df,data3)
}

save(obs_df,file='./data processed/obs_depth_density.RData')

obs_df$species<-factor(obs_df$species, levels = spp1)

# Kernel density estimation for depth distributions
ggplot(obs_df, aes(x = depth_m, fill = year_type)) +
  geom_density(alpha = 0.5) +
  #facet_wrap(~ Species) +
  labs(x = 'Depth (m)', y = 'Density', title = 'Depth Density Distribution in Warm vs. Cold Years') +
  scale_fill_manual(values = c('warm' = 'red', 'cold' = 'blue'), name = 'Year Category') +
  theme_minimal()+
  #scale_x_continuous(limits = c(0,600))+
  theme(strip.text = element_text(size=12),plot.title = element_text(hjust=0.5))+
  facet_wrap(~species)




#Get observed data temp, depth and biomass for 

#####################
# PREDICTED DATA
#####################

spp<-list.files('./data processed/species/')
spp1<-spp[c(1:24)]

dens_all<-data.frame(matrix(NA,nrow=0,ncol=13))
colnames(dens_all)<-c("Lat","Lon","Area_in_survey_km2" ,"Stratum","region","DepthGEBCO","depth_m","Temp","Year","dens.Site",         
                    "dens.Time","dens.value",'species')

ind_all<-data.frame(matrix(NA,nrow=0,ncol=5))
colnames(ind_all)<-c("year","EBSshelf", "NBS","EBSslope" ,"species" )

# Loop over species
for (s in sp_shelfslope[-6]) { #[c(1,2,4,7)]
  
  #s<-sp_shelfslope[4]
  
  cat(paste('#################### ',s," - ",'\n'))
  
  # Load shelf model
  load(paste0('./shelf EBS NBS VAST/', s, '/fit.RData'))
  fit_shelf <- fit
  bio_shelf <- fit_shelf$Report$Index_gctl[,,,1]
  x<-fit_shelf$Report$D_gct[,1,] #dens
  dens<-reshape2::melt(drop_units(x))
  ind_NBS<-fit_shelf$Report$Index_ctl[,,2] #index NBS
  ind_EBSshelf<-fit_shelf$Report$Index_ctl[,,3] #index EBSshelf
  
  if (s=='Atheresthes evermanni') {
    dens_shelf<-
      data.frame(
        grid.ebs_year[which(grid.ebs_year$region %in% c('EBSshelf','NBS') & grid.ebs_year$Year %in% c(1991:2022)),],
        dens=dens)
    
    ind_shelf<-
      data.frame(year=1991:2022,
                 EBSshelf=drop_units(ind_EBSshelf),
                 NBS=drop_units(ind_NBS))
    
  } else {
    dens_shelf<-
      data.frame(
        grid.ebs_year[which(grid.ebs_year$region %in% c('EBSshelf','NBS') & grid.ebs_year$Year %in% c(1982:2022)),],
        dens=dens)
    
    ind_shelf<-
      data.frame(year=1982:2022,
                 EBSshelf=drop_units(ind_EBSshelf),
                 NBS=drop_units(ind_NBS))
    
  }

  
  # Load slope model
  load(paste0('./slope EBS VAST/', s, '/fit_st.RData'))
  fit_slope <- fit
  bio_slope <- fit_slope$Report$Index_gctl[,,,1]
  x<-fit_slope$Report$D_gct[,1,] #dens
  dens<-reshape2::melt(drop_units(x))
  ind_slope<-fit_slope$Report$Index_ctl[,,1] #index slope
  
  dens_slope<-
    data.frame(
      grid.ebs_year[which(grid.ebs_year$region %in% c('EBSslope') & grid.ebs_year$Year %in% c(2002:2016)),],
      dens=dens)
  
  dens_df<-rbind(dens_shelf,dens_slope)
  
  dens_df$species<-s
  
  dens_all<-rbind(dens_all,dens_df)
  
  if (s=='Atheresthes evermanni') {
    
    ind_df<-
      data.frame(rbind(data.frame('year'=1982:1990,
                                  'EBSshelf'=c(rep(NA,times=length(1982:1990))),
                                  'NBS'=c(rep(NA,times=length(1982:1990)))),
                       ind_shelf),
                 EBSslope=c(rep(NA,times=20),drop_units(ind_slope),rep(NA,times=6)))

  } else {
    
    ind_df<-
      data.frame(ind_shelf,
                 EBSslope=c(rep(NA,times=20),drop_units(ind_slope),rep(NA,times=6)))
  }
  
  ind_df$species<-s
  
  ind_all<-rbind(ind_all,ind_df)
  
  
}
  




names(dens_all)<-c("Lat","Lon","Area_in_survey_km2" ,"Stratum","region","DepthGEBCO","depth_m","Temp","Year","Site",         
                   "Time","dens","species")

save(dens_all,file='./data processed/dens_shelfslope.RData')
save(ind_all,file='./data processed/ind_shelfslope.RData')

#cold and warm years
wyrs<-c(2016,2018,2022)
cyrs<-c(1999,2008,2009,2010,2012)

#dens_all1<-subset(dens_all,dens.value!=0 ) #& region!='NBS'
dens_all1<-dens_all[which(dens_all$dens>=0.00001 & dens_all$region!='NBS' & dens_all$Year %in% c(cyrs,wyrs)),] #& region!='NBS' 
summary(dens_all1$dens)
dim(dens_all1)
dim(dens_all)


dens_all1$year_type<-ifelse(dens_all1$Year %in% wyrs, "warm",'cold')

# Kernel density estimation for depth distributions
ggplot(dens_all1, aes(x = depth_m, fill = year_type)) +
  geom_density(alpha = 0.5) +
  #facet_wrap(~ Species) +
  labs(x = 'Depth (m)', y = 'Density', title = 'Depth Density Distribution in Warm vs. Cold Years') +
  scale_fill_manual(values = c('warm' = 'red', 'cold' = 'blue'), name = 'Year Category') +
  theme_minimal()+
  scale_x_continuous(limits = c(0,600))+
  theme(strip.text = element_text(size=12),plot.title = element_text(hjust=0.5))+
  facet_wrap(~species)



###################

ind_all1<-reshape2::melt(ind_all,id.vars=c('year','species'))

#remove 'all' and 'EBS+NBS' region
ind_all2<-subset(ind_all1,variable %in% c('NBS','EBSshelf'))

ind_all3<-
  ind_all2 %>% group_by(species, year) %>%
  mutate(frac = value / sum(value))

#ind_all3$year<-as.factor(ind_all3$year)

# Your original plot code with correction for coloring axis text
#p<-
ggplot() +
  annotate("rect", xmin = 2002, xmax = 2005.5, ymin = -Inf, ymax = Inf, alpha = 0.5, fill = 'red') +
  annotate("rect", xmin = 2013.5, xmax = 2022, ymin = -Inf, ymax = Inf, alpha = 0.5, fill = 'red') +
  annotate("rect", xmin = 2005.5, xmax = 2013.5, ymin = -Inf, ymax = Inf, alpha = 0.5, fill = 'blue') +
  geom_bar(data = ind_all3, aes(x = year, y = frac, fill = variable), stat = 'identity') +
  scale_fill_manual(values = c('EBSshelf' = '#046407', 'NBS' = '#B4AF46'),
                    breaks = c('EBSshelf', 'NBS'),
                    labels = c('EBSshelf', 'NBS'),name='region') +
  scale_x_continuous(breaks = c(1985, 1990, 1995, 2000, 2005, 2010, 2015, 2020)) +
  labs(y = 'abundance proportion') +
  theme_bw() +
  theme(text = element_text(size = 14)) +

  facet_wrap(~ species, scales = 'free_x', nrow = 2)

  #theme(axis.text.x = element_text(angle = 45,vjust=0.5,hjust=0.5,color = ifelse(levels(ind_all3$year) %in% c(2003:2005,2015:2016), "red", 
   #                                               ifelse(levels(ind_all3$year) %in% c(2009,2008,2012), "blue",'black'))))






