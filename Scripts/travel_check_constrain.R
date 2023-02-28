library(dplyr)
library(lubridate)
library(geosphere)
library(TSP) #Traveling Salesperson Problem

#read data
x<-readRDS('C:/Users/Daniel.Vilas/Work/Adapting Monitoring to a Changing Seascape/data raw/afsc_haul_raw_2023_2_21.rds')

x$year<-year(as.POSIXlt(x$date, format="%d/%m/%Y"))

x %>% count(survey_name,year)


xx<-subset(x,survey_name=='Eastern Bering Sea Crab/Groundfish Bottom Trawl Survey')

years<-unique(xx$year)



for (y in years) {
  
 y<-2021 
  
xxx<-subset(xx,year==y)


df<-data.frame('ID'=xxx$hauljoin,'Lat'=xxx$lat_start,'Lon'=xxx$lon_start)

ggplot(data=df,aes(x=Lon,y=Lat))+
  geom_point()+
  geom_text(aes(label=ID),vjust=-0.5)

m<-distm(df[3:2],df[3:2])

# bind the distance matrix to the dataframe
df = cbind.data.frame(df, m)
 
#solve for shortest path 
tsp = TSP(m)
tour = solve_TSP(tsp)

#tour length 
tour_length(tour)

#permutation vector for shortest tour
perm_vec = as.integer(tour)
 
#re-arrange cities data frame in order of city tour
df <- df[perm_vec,]
df$order<-c(1:nrow(df))
#label ID and n order
df$label<-paste0(df$ID,'n',df$order)

#plot
ggplot()+
  #geom_text(data=df,aes(x=Lon,y=Lat,label=label))+
  geom_path(data=df,aes(x=Lon,y=Lat))
  

}














      