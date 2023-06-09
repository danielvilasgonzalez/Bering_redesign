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
                     n_strata=c(15),stringsAsFactors = FALSE) #c(5,10,15)

#add scenario number
samp_df$samp_scn<-paste0(paste0('scn',1:nrow(samp_df)))
samp_df<-rbind(samp_df,c('baseline','current',520,15,'scnbase'),
                         c('baseline w/o corner','current',494,15,'scnbase_bis'))
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
  
  #print scenario to check progress
  cat(paste(" ##########################################################################\n",
            " #############   sampling scheme", samp_df[s,'samp_scn'],  "  #############\n",
            " ##########################################################################\n"))

  if (grepl('base',samp_df[s,'samp_scn'])) { #if it contains base is baseline scenario
    
    strata<-baseline_strata$cell_strata[,c('cell','Lat','Lon','Stratum')]
    names(strata)[4]<-c('Strata')
    D8<-strata
    
    #allocations<-
      if(samp_df[s,'samp_scn']=='scnbase'){
        allocations<-data.frame('Strata'=baseline_strata$n_samples$stratum,
                                'n_samples'=baseline_strata$n_samples$scnbase)
        }else{
        allocations<-data.frame('Strata'=baseline_strata$n_samples$stratum,
                                'n_samples'=baseline_strata$n_samples$scnbase_bis)
        }

  } else {
  
  #load results_optimization
  load(file=paste0('./output/species/',sp,'/optimization data/optimization_results_',samp_df[s,'samp_scn'],'.RData')) #result_list
  strata<-result_list$solution$indices
  colnames(strata)<-c('cell','Strata')
  allocations<-data.frame('Strata'=1:length(result_list$sample_allocations),
                          'n_samples'=result_list$sample_allocations)
  
  #add a strata value to each cell
  D8<-merge(D6,strata,by='cell',all.x=TRUE)
  D8<-D8[,c("cell","Lat","Lon","Strata")]
  D8$Strata<-as.numeric(D8$Strata)
  D8$Strata<-ifelse(is.na(D8$Strata),999,D8$Strata)
  
  }
  
  #df to spatialpoint df
  coordinates(D8) <- ~ Lon + Lat
  crs(D8)<-c(crs="+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
  
  #reproject coordinates for plotting purposes
  D8_1<-spTransform(D8,'+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs')
  D8_2<-merge(as.data.frame(D8_1),x5,by='cell')
  D8_2<-D8_2[,c('cell','Strata','Lon.x','Lat.x','col','row')]
  names(D8_2)[c(2,3,4)]<-c('strata','Lon','Lat')
  
  #to store samples
  all_points<-array(NA,
                    dim = list(sum(allocations$n_samples),4,length(1982:2027),3),
                    dimnames = list(c(1:sum(allocations$n_samples)),c('Lon','Lat','cell','strata'),1:length(1982:2027),c('current','buffer','random')))

  for (y in 1:length(1982:2027)) {
    
      #y<-1
      
      #print scenario to check progress
      cat(paste(" #############   year", y, 'of',length(1982:2027),  "  #############\n"))
      
      #to store points
      dfcurrent<-data.frame(matrix(NA,nrow=0,ncol=4))
      colnames(dfcurrent)<-c('Lon','Lat','cell','strata')
      dfbuffer<-dfrandom<-dfcurrent
      
      #for while purposes
      flag<-TRUE
      
      #random sample for each strata wit distance constrains
      for(n_istrata in 1:nrow(allocations)) {
        
        #n_istrata<-10
        
        istrata<-allocations[n_istrata,'Strata']
      
        #subset cells for strata
        df<-subset(as.data.frame(D8_2),strata==istrata)
        df1<-subset(as.data.frame(D8_2),strata==istrata & cell %in% baseline_strata$locations$cell)
        n_i<-allocations[n_istrata,'n_samples']
        
        #while error, keep running
        while(flag){
           
           #keep running if error
           tryCatch({
            
             #print loop state
             cat(paste(" -- strata", n_istrata, 'of',nrow(allocations),  " -- CURRENT  --\n"))
             
        ##############################
        # CURRENT STRATIFIED APPROACH
        ##############################
        
        if(samp_df[s,'samp_scn']=='scnbase'){
              pointsc<-baseline_strata$locations[which(baseline_strata$locations$stratum==istrata),c('longitude','latitude','cell','stratum')]
              names(pointsc)<-c('Lon' ,'Lat', 'cell','strata')
        }else if (samp_df[s,'samp_scn']=='scnbase_bis'){
              pointsc<-baseline_strata$locations[which(baseline_strata$locations$corner==FALSE & baseline_strata$locations$stratum==istrata),c('longitude','latitude','cell','stratum')]
              names(pointsc)<-c('Lon' ,'Lat',  'cell','strata')
        } else {
          
        #if more required samples than available
        if (nrow(df1)<n_i) {
          
          #vectors to store results
          dropcell<-c()
          selcell<-c()
          
          #duplicate df strata
          dff<-df
          
          #df removing available samples from current design
          dff<-subset(dff, !(cell %in% df1$cell))
          
          #cells to complete the required cells
          ii<-n_i-nrow(df1)
          
          #loop over the required samples using buffer
          for (iii in rep(1,times=ii)) {
            
              #ii<-1
              
              #get random cell
              cell_i<-sample(dff$cell,iii)
              #get row
              #row<-dff[which(dff$cell==cell_i),c('row','col')][1,1] 
              row<-dff[which(dff$cell %in% c(cell_i,df1$cell)),c('row','col')][1,1]
              #get col
              #col<-dff[which(dff$cell==cell_i),c('row','col')][1,2] 
              col<-dff[which(dff$cell %in% c(cell_i,df1$cell)),c('row','col')][1,2]
              #get adjacent cells
              adj_cells<-expand.grid(row=c(row,row+(1:100),row-(1:100)),
                                     col=c(col,col+(1:100),col-(1:100)))
              adj_cells1<-merge(adj_cells,df,by=c('row', 'col'))[,'cell']
              #remove cells from available for sampling
              dropcell_i<-c(cell_i,adj_cells1)

              #if no more samples available to sample stop
              if (nrow(subset(dff, !(cell %in% dropcell_i)))==0) {flag<-TRUE; stop("-- no samples to select on current approach",call. = FALSE)}else {
                dff<-subset(dff, !(cell %in% dropcell_i))
                selcell<-c(selcell,cell_i)
                dropcell<-c(dropcell,dropcell_i)}
          }
          
          #subset points if equal to required number of samples, if not ERROR and rerun iteration
          if (nrow(subset(df,cell %in% c(selcell,df1$cell))) != n_i) {flag<-TRUE; stop("-- different number of points on current approach",call. = FALSE)}else{
            pointsc<-subset(df,cell %in% c(selcell,df1$cell))[,c('Lon','Lat','cell','strata')]
            names(pointsc)<-c('Lon','Lat','cell','strata')}
        
        #else there are enough samples to get from the current sampling design    
        } else {
          ss<-sample(1:nrow(df1),size = n_i,replace = FALSE)
          pointsc<-df1[ss,c('Lon','Lat','cell','strata')]}}
             
        ##############################
        # STRATIFIED + BUFFER APPROACH
        ##############################
        
        cat(paste(" -- strata", n_istrata, 'of',nrow(allocations),  " -- BUFFER  --\n"))

        #replicate df of strata     
        dff<-df
        
        #create vectors to store cells
        dropcell<-c()
        selcell<-c()
        
        #loop over required samples
        for (iii in rep(1,times=n_i)) {
          
          #iii<-1
          
          #random sample
          cell_i<-sample(dff$cell,iii)
          #get row of selected sample
          row<-dff[which(dff$cell==cell_i),c('row','col')][1,1]
          #get col of selected sample
          col<-dff[which(dff$cell==cell_i),c('row','col')][1,2]
          #get adjacent cells of selected sample
          adj_cells<-expand.grid(row=c(row,row+(1:100),row-(1:100)),
                                 col=c(col,col+(1:100),col-(1:100)))
          adj_cells1<-merge(adj_cells,df,by=c('row', 'col'))[,'cell']
          #drop cells is equal to selected cell and adjacent cells
          dropcell_i<-c(cell_i,adj_cells1)
          
          #store cells if there are still samples on the filtered df
          if (nrow(subset(dff, !(cell %in% dropcell_i)))==0) {flag<-TRUE;stop("-- no samples to select on buffer approach",call. = FALSE)} else {
              dff<-subset(dff, !(cell %in% dropcell_i))
              selcell<-c(selcell,cell_i)
              dropcell<-c(dropcell,dropcell_i)}

        }
        
        #subset points if equal to required number of samples, if not ERROR and rerun iteration
        if (nrow(subset(df,cell %in% selcell)) != n_i ) {flag<-TRUE;stop("-- different number of points on buffer approach",call. = FALSE)}else{
          pointsb<-subset(df,cell %in% selcell)[,c('Lon','Lat','cell','strata')]
          names(pointsb)<-c('Lon','Lat','cell','strata')}
        
        ##############################
        # STRATIFIED RANDOM APPROACH
        ##############################
        
        cat(paste(" -- strata", n_istrata, 'of',nrow(allocations),  " -- RANDOM  --\n"))
        
        #duplicate df strata
        dff<-df
        
        #random selection of samples
        selcell<-sample(dff$cell,n_i)
        pointsr<-subset(df,cell %in% selcell)[,c('Lon','Lat','cell','strata')]
        names(pointsr)<-c('Lon','Lat','cell','strata')
        
        #append data if buffer and current df have equal to required number of samples
        if (nrow(pointsc) == n_i & nrow(pointsb) == n_i & nrow(pointsr) == n_i) {
          dfcurrent<-rbind(dfcurrent,pointsc)
          dfbuffer<-rbind(dfbuffer,pointsb)
          dfrandom<-rbind(dfrandom,pointsr)
          flag<-FALSE}
        },
        
        #error message
        error=function(e) {
          message(paste0("ERROR RETRYING",':\n'),e)})
          if (!flag) next
        }
        
        #to start while again
        flag<-TRUE
      }
      
      #append points
      all_points[,,y,'current']<-unlist(dfcurrent)
      all_points[,,y,'buffer']<-unlist(dfbuffer)
      all_points[,,y,'random']<-unlist(dfrandom)
    }
  
#save sample selection
save(all_points,file=paste0('./output/species/',sp,'/optimization data/samples_optimization_',samp_df[s,'samp_scn'],'_dynamic.RData'))
}
