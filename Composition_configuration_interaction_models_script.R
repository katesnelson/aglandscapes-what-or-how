library(pacman)

p_load(tidyverse, INLA, sf, spdep, rgeos, knitr, ggpubr,doParallel, foreach,tidyselect, viridis)


wd<-getwd()

source(paste0(wd,'/Scripts/SA_functions.R'))



########################################################
### Run composition-configuration interaction models ###
########################################################

#Build corn data    
  data.quants(file="corn_panel", savename=paste0("corn_","quants"), crop="corn", 
              projection="+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m no_defs") #102003 epsg code for USA albers equal area conic


  #SHDIxED model
    f1<-build.formula.facint1("LSM_SHDI_ALL", "corn", controls = "All", priors = c("default")) 
    
    m1<- run.model(f1, data=paste0("corn_","quants"), crop="corn", priors ="default", name=paste0("corn","_","SHDIxEDquants"))
    
  
  #SHDIxMN_AREA model
    f1<-build.formula.facint2("LSM_SHDI_ALL", "corn", controls = "All", priors = c("default")) 
    
    m1<- run.model(f1, data=paste0("corn_","quants"), crop="corn", priors ="default", name=paste0("corn","_","SHDIxMNAREAquants"))
  
  
  #RICHxED model
    f1<-build.formula.facint1("LSM_RICH_ALL", "corn", controls = "All", priors = c("default")) 
    
    m1<- run.model(f1, data=paste0("corn_","quants"), crop="corn", priors ="default", name=paste0("corn","_","RICHxEDquants"))
  
  
  #RICHxMN_AREA model
    f1<-build.formula.facint2("LSM_RICH_ALL", "corn", controls = "All", priors = c("default")) 
    
    m1<- run.model(f1, data=paste0("corn_","quants"), crop="corn", priors ="default", name=paste0("corn","_","RICHxMNAREAquants"))


#Build soy data
  data.quants(file="soy_panel", savename=paste0("soy_","quants"), crop="soy", 
              projection="+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m no_defs") #102003 epsg code for USA albers equal area conic


  #SHDIxED model
    f1<-build.formula.facint1("LSM_SHDI_ALL", "soy", controls = "All", priors = c("default")) 
    
    m1<- run.model(f1, data=paste0("soy_","quants"), crop="soy", priors ="default", name=paste0("soy","_","SHDIxEDquants"))
  
  
  #SHDIxMN_AREA model
    f1<-build.formula.facint2("LSM_SHDI_ALL", "soy", controls = "All", priors = c("default")) 
    
    m1<- run.model(f1, data=paste0("soy_","quants"), crop="soy", priors ="default", name=paste0("soy","_","SHDIxMNAREAquants"))
  
  
  #RICHxED model
    f1<-build.formula.facint1("LSM_RICH_ALL", "soy", controls = "All", priors = c("default")) 
    
    m1<- run.model(f1, data=paste0("soy_","quants"), crop="soy", priors ="default", name=paste0("soy","_","RICHxEDquants"))
  
  
  #RICHxMN_AREA model
    f1<-build.formula.facint2("LSM_RICH_ALL", "soy", controls = "All", priors = c("default")) 
    
    m1<- run.model(f1, data=paste0("soy_","quants"), crop="soy", priors ="default", name=paste0("soy","_","RICHxMNAREAquants"))



#Build wheat data
  data.quants(file="wwheat_panel", savename=paste0("wwheat_","quants"), crop="wwheat", 
              projection="+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m no_defs") #102003 epsg code for USA albers equal area conic


  #SHDIxED model
    f1<-build.formula.facint1("LSM_SHDI_ALL", "wwheat", controls = "All", priors = c("default")) 
    
    m1<- run.model(f1, data=paste0("wwheat_","quants"), crop="wwheat", priors ="default", name=paste0("wwheat","_","SHDIxEDquants"))
  
  
  #SHDIxMN_AREA model
    f1<-build.formula.facint2("LSM_SHDI_ALL", "wwheat", controls = "All", priors = c("default")) 
    
    m1<- run.model(f1, data=paste0("wwheat_","quants"), crop="wwheat", priors ="default", name=paste0("wwheat","_","SHDIxMNAREAquants"))
  
  
  #RICHxED model
    f1<-build.formula.facint1("LSM_RICH_ALL", "wwheat", controls = "All", priors = c("default")) 
    
    m1<- run.model(f1, data=paste0("wwheat_","quants"), crop="wwheat", priors ="default", name=paste0("wwheat","_","RICHxEDquants"))
  
  
  #RICHxMN_AREA model
    f1<-build.formula.facint2("LSM_RICH_ALL", "wwheat", controls = "All", priors = c("default")) 
    
    m1<- run.model(f1, data=paste0("wwheat_","quants"), crop="wwheat", priors ="default", name=paste0("wwheat","_","RICHxMNAREAquants"))



