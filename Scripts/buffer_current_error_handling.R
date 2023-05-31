#libraries from cran to call or install/load
pack_cran<-c("splines",'SamplingStrata','wesanderson','dplyr','sp',
             'sf','maptools','parallel','rasterVis','rgeos','scales',
             'rnaturalearth','grid','ggplot2','spatstat')

#install pacman to use p_load function - call library and if not installed, then install
if (!('pacman' %in% installed.packages())) {
  install.packages("pacman")}

#install akgfmaps to extract shapefile of Alaska
if (!('akgfmaps' %in% installed.packages())) {
  devtools::install_github('afsc-gap-products/akgfmaps')};library(akgfmaps)

#install coldpool to extract SBT for the EBS
if (!('VAST' %in% installed.packages())) {
  devtools::install_github("james-thorson/VAST@main", INSTALL_opts="--no-staged-install")};library(VAST)

#load/install packages
pacman::p_load(pack_cran,character.only = TRUE)

#call function distance points buffer
#source('C:/Users/Daniel.Vilas/Work/GitHub/Bering_redesign/Scripts/genRandomPnts.R')

#setwd
out_dir<-'C:/Users/Daniel.Vilas/Work/Adapting Monitoring to a Changing Seascape/'
setwd(out_dir)

###################################
# LOAD GRID EBS (remember to keep the same order as in fit_model if multiple grids)
###################################

#load grid of NBS and EBS
load('./extrapolation grids/northern_bering_sea_grid.rda')
load('./extrapolation grids/eastern_bering_sea_grid.rda')
grid<-as.data.frame(rbind(data.frame(northern_bering_sea_grid,region='NBS'),data.frame(eastern_bering_sea_grid,region='EBS')))
grid$cell<-1:nrow(grid)
#add col and row number
x1<-grid[,c('Lon','Lat','cell')]
names(x1)<-c('x','y','z')
coordinates(x1)=~x + y
crs(x1)<-c(crs="+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
x2<-spTransform(x1,'+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs')
x3<-data.frame(x2)
x3$x<-as.integer(x3$x)
x3$y<-as.integer(x3$y)
lon<-sort(unique(x3$x),decreasing = FALSE) #1556
lat<-sort(unique(x3$y),decreasing = TRUE) #1507
lons<-data.frame(x=lon,col=1:length(lon))
lats<-data.frame(y=lat,row=1:length(lat))
x4<-merge(x3,lons,by='x',all.x=TRUE)
x5<-merge(x4,lats,by='y',all.x=TRUE)
colnames(x5)<-c('Lat','Lon','cell','optional','col','row')
grid<-x5[,c('Lat','Lon','cell','col','row')]

###################################
# SCENARIOS
###################################

#sampling scenarios
samp_df<-expand.grid(strat_var=c('Depth_varTemp','varTemp','Depth'),
                     target_var=c('sumDensity'), #,'sqsumDensity'
                     n_samples=c(520), #c(300,500) 520 (EBS+NBS+CRAB);26 (CRAB); 350 (EBS-CRAB); 494 (NBS-CRAB)
                     n_strata=c(15)) #c(5,10,15)

#add scenario number
samp_df$samp_scn<-paste0(paste0('scn',1:nrow(samp_df)))

load('./output/baseline_strata.RData') #baseline_strata
baseline_strata$locations[grep('GF|HG|JI|IH|ON|QP|PO',baseline_strata$locations$stationid),]
baseline_strata$locations$corner<-ifelse(grepl('GF|HG|JI|IH|ON|QP|PO',baseline_strata$locations$stationid),'TRUE','FALSE')

###################################
# LOOP
###################################

sp<-'Gadus macrocephalus'

#load optimization data
load(paste0('./output/species/',sp,'/optimization data/optimization_static_data.RData')) #D6
#load(paste0('./output/species/',sp,'/projection_data.RData')) #temp_dens_vals

#load fit OM
load(paste0('./shelf EBS NBS VAST/',sp,'/fit.RData'))

#loop through sampling scenarios
for (s in 1:nrow(samp_df)) {

s<-3
    
#load results_optimization
load(file=paste0('./output/species/',sp,'/optimization data/optimization_results_',samp_df[s,'samp_scn'],'.RData')) #result_list
strata<-result_list$solution$indices
colnames(strata)<-c('cell','Strata')
allocations<-result_list$sample_allocations


#add a strata value to each cell
D8<-merge(D6,strata,by='cell',all.x=TRUE)
D8<-D8[,c("cell","Lat","Lon","Strata")]
D8$Strata<-as.numeric(D8$Strata)
D8$Strata<-ifelse(is.na(D8$Strata),999,D8$Strata)

#df to spatialpoint df
coordinates(D8) <- ~ Lon + Lat
crs(D8)<-c(crs="+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")

#reproject coordinates for plotting purposes
D8_1<-spTransform(D8,'+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs')
D8_2<-merge(as.data.frame(D8_1),x5,by='cell')
D8_2<-D8_2[,c('cell','Strata','Lon.x','Lat.x','col','row')]
names(D8_2)[c(2,3,4)]<-c('strata','Lon','Lat')
  
#iterations for samples
n_iter<-100

#to store samples
all_points<-array(NA,dim = list(sum(allocations),4,n_iter,2),
                  dimnames = list(c(1:sum(allocations)),c('Lon','Lat','cell','strata'),c(1:n_iter),c('current','buffer')))

flag<-TRUE

#loop over iterations
for (iter in 1:n_iter) {
  
  while(flag){
    
    tryCatch({
  
  #print scenario to check progress
  cat(paste(" #############   iter", iter, 'of',n_iter,  "  #############\n"))
 
  #iter<-1
  
  #to store points
  dfcurrent<-data.frame(matrix(NA,nrow=0,ncol=4))
  colnames(dfcurrent)<-c('Lon','Lat','cell','strata')
  dfbuffer<-dfcurrent
  
  #random sample for each strata wit distance constrains
  for(istrata in 1:length(allocations)) {
    #istrata<-1
    
    ##############################
    # CURRENT APPROACH
    ##############################
    
    #subset cells for strata
    df<-subset(as.data.frame(D8_2),strata==istrata)
    df1<-subset(as.data.frame(D8_2),strata==istrata & cell %in% baseline_strata$locations$cell)
    
    n_i<-allocations[istrata]
    
    if (nrow(df1)<n_i) {
      
      dropcell<-c()
      selcell<-c()
      
      dff<-df
      
      ii<-n_i-nrow(df1)
      

      
      for (ii in rep(1,times=n_i)) {
        
        
            #ii<-1
            cell_i<-sample(dff$cell,ii)
            row<-dff[which(dff$cell==cell_i),c('row','col')][1,1]
            col<-dff[which(dff$cell==cell_i),c('row','col')][1,2]
            adj_cells<-expand.grid(row=c(row,row+(1:200),row-(1:200)),
                                   col=c(col,col+(1:200),col-(1:200)))
            adj_cells1<-merge(adj_cells,df,by=c('row', 'col'))[,'cell']
            # if (istrata == 2 ) {
            #   x<-merge(adj_cells,df,by=c('row', 'col'))
            # }
            
            dropcell_i<-c(cell_i,adj_cells1)
            selcell<-c(selcell,cell_i)
            dropcell<-c(dropcell,dropcell_i)
            
            dff<-subset(dff, !(cell %in% dropcell))
                    
        
        
      }
      
      points<-subset(df,cell %in% selcell)[,c('Lon','Lat','cell','strata')]
      names(points)<-c('Lon','Lat','cell','strata')
      
      #append data
      dfcurrent<-rbind(dfcurrent,points)
      
    } else {
      
      ss<-sample(1:nrow(df1),size = allocations[istrata],replace = FALSE)
      dfcurrent<-rbind(dfcurrent,df1[ss,c('Lon','Lat','cell','strata')])
      
    }
    
    ##############################
    # BUFFER APPROACH
    ##############################
    
    dff<-df
    
    dropcell<-c()
    selcell<-c()
    
    
    for (ii in rep(1,times=n_i)) {
      
     
          
          #ii<-1
          cell_i<-sample(dff$cell,ii)
          row<-dff[which(dff$cell==cell_i),c('row','col')][1,1]
          col<-dff[which(dff$cell==cell_i),c('row','col')][1,2]
          adj_cells<-expand.grid(row=c(row,row+(1:200),row-(1:200)),
                                 col=c(col,col+(1:200),col-(1:200)))
          adj_cells1<-merge(adj_cells,df,by=c('row', 'col'))[,'cell']
          # if (istrata == 2 ) {
          #   x<-merge(adj_cells,df,by=c('row', 'col'))
          # }
          
          dropcell_i<-c(cell_i,adj_cells1)
          selcell<-c(selcell,cell_i)
          dropcell<-c(dropcell,dropcell_i)
          
          dff<-subset(dff, !(cell %in% dropcell))
      
    }
    
    # if (length(selcell)!=i) {
    #   #print scenario to check progress
    #   cat(paste(" #############  not enough samples  #############\n"))}
    
    points<-subset(df,cell %in% selcell)[,c('Lon','Lat','cell','strata')]
    names(points)<-c('Lon','Lat','cell','strata')
    
    #append data
    dfbuffer<-rbind(dfbuffer,points)
    
  }
  
  #append points
  all_points[,,iter,'current']<-unlist(dfcurrent)
  all_points[,,iter,'buffer']<-unlist(dfbuffer)
  
  flag<-FALSE},
  error=function(e) {
    
    
    message(paste0("Daniel tienes un problema",':\n'),e)})
    
    if (!flag) next
    
  }
  
  flag<-TRUE
  
}


#save sample selection
save(all_points,file=paste0('./output/species/',sp,'/optimization data/samples_optimization_',samp_df[s,'samp_scn'],'.RData'))


}
