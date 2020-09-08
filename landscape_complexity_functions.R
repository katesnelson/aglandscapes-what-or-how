
#######################################################
###function to spatial join by largest spatial overlap
######################################################
area.overlay<-function(file1, file2, projection,id1,id2, savename){

  d1<-readRDS(paste0(wd,"/data/",file1,".RDS"))
  d1<-st_as_sf(d1)
  if (st_crs(d1)[[2]]!= projection){
    d1<-st_transform(d1, projection)
  }
  d1<-d1[,c(paste0(id1),"geometry")]
  
  
  d2<-readRDS(paste0(wd,"/data/",file2,".rds"))
  d2<-d2 %>% st_as_sf(.)%>% st_transform(.,projection) #use an equal-area projection for U.S.
  d2$AERCODE<-seq(1:length(d2$US_L3CODE)) #ONLY for Ag-diversity project --> build a unique AERCODE number for each AER polygon (so discontinuous areas with same AER class are considered to be different)
  d2<-d2[,c(paste0(id2),"geometry")]
  
  #join Counties to Agro-Eco-Regions by assigning Counties to the AER with which they have the largest spatial overlap
  d1_d2<- st_intersection(d1,d2) %>% mutate(int_area = as.numeric(st_area(.))) %>% arrange(., desc(int_area)) #sort from large to small intersection area
  d3<-d1_d2[duplicated(d1_d2[,paste0(id1)][[1]])==FALSE, ]
  saveRDS(d3, paste0(wd,"/Data/",savename,".RDS"))

}


#######################################
###function to prep data for modeling
######################################
prep.data<-function(file, savename, crop, projection){
 
  dat<-readRDS(paste0(wd,"/Data/",file,".RDS")) #read in the data
  dat<-dat %>% mutate(., Yr = as.numeric(as.character(YEAR))) %>% mutate(., Yr= Yr - min(Yr)) #this will set the first year as the base year (intercept will be based on this year)
  dat$GEOID<-as.character(dat$GEOID)
  
  cnty<-readRDS(paste0(wd,"/Data/county.RDS"))
  cnty<-st_as_sf(cnty) %>% dplyr::select(.,c(5,11)) %>% st_transform(., projection)
  
  cnty_AER<-readRDS(paste0(wd,"/Data/cnty_AER_SA.RDS")) %>% st_transform(., projection)
  
  #join to add county and AER ids to dataset
  d<-left_join(cnty,dat, by="GEOID")
  d<-left_join(d,st_set_geometry(cnty_AER[,c(1,2)], NULL), by='GEOID')
  
  
  #prep for building adjacency matrix
  d_sub<-d[!is.na(d$YIELD),]
  names<-unique(d_sub$GEOID)

  #clip county file to crop producing areas only and setup indexing by CNTY
  cnty<-cnty[cnty$GEOID %in% names,] %>% arrange(., GEOID) %>% mutate (.,CNTY=seq(1,nrow(.),1)) #order dataset by county, build group index that corresponds to adjacency matrix

  #Build the relational matrix for areas of interest
  neighbors<-as(cnty,"Spatial") %>% poly2nb(., queen=F)#convert county sf to spatial polygons for inla relation matrix & create neighbors list from polygon object (neighbors share one or points at boundary)
  H.adj <- nb2mat(neighbors, style ="B", zero.policy=TRUE ) #convert to a sparse matrix to reduce memory (neighbor list to binary coded neighbor weights matrix)
  H.adj <-as(H.adj, "dgTMatrix") #convert to a sparse style matrix
  saveRDS(H.adj, paste0(wd,"/Data/H.adj.",crop,".rds"))
  
  #plot spatial object
  #image(inla.graph2matrix(H.adj), xlab="", ylab="")
  
  #add indexes and set up vars for log-linear model with rw (by binning) and quadratic (by squaring) terms
  d_sub<-left_join(d_sub,st_set_geometry(cnty[,c("GEOID","CNTY")], NULL),by="GEOID")#join new index to crop dataset
  d_sub<-d_sub %>% arrange(.,CNTY,YEAR) #Ordering dataset by county then time for modeling 
  
  d_sub$AERCODE.id<-d_sub$AERCODE
  d_sub$YEAR.id<-d_sub$Yr  
  d_sub$YEAR.id2<-(d_sub$Yr )^2
  d_sub$AERCODE.id2<-d_sub$AERCODE
  
  #Bin data for rw1 evaluation
  d_sub<- d_sub %>% mutate_at(ivars, ~inla.group(., n = 20, method = "quantile")) %>%   
    mutate_at(c("TP","GDD","SDD", "SOIL"), ~inla.group(., n = 20, method = "quantile")) 
  
  d_sub$PERC_IRR[is.na(d_sub$PERC_IRR)]<-0
  
  #Yield transformed to log(Yield)
  d_sub$YIELD<-log(d_sub$YIELD)
  
  
  saveRDS(d_sub,paste0(wd,"/Data/",savename,".rds"))
  
 }


#####################################################
###function to prep standardized data for modeling (same as prep.data but indep vars are standardized)
###################################################
prep.data.std<-function(file, savename, crop, projection){
  
  dat<-readRDS(paste0(wd,"/Data/",file,".RDS")) #read in the data
  dat<-dat %>% mutate(., Yr = as.numeric(as.character(YEAR))) %>% mutate(., Yr= Yr - min(Yr)) 
  dat$GEOID<-as.character(dat$GEOID)
  
  cnty<-readRDS(paste0(wd,"/Data/county.RDS"))
  cnty<-st_as_sf(cnty) %>% dplyr::select(.,c(5,11)) %>% st_transform(., projection)
  
  cnty_AER<-readRDS(paste0(wd,"/Data/cnty_AER_SA.RDS"))%>% st_transform(., projection)
  
  #join to add county and AER ids to dataset
  d<-left_join(cnty,dat, by="GEOID")
  d<-left_join(d,st_set_geometry(cnty_AER[,c(1,2)], NULL), by='GEOID')
  
  # t<-d %>% mutate (AERCODE.y=as.factor(AERCODE.y))
  # plot(t[t$Yr==0,"AERCODE.y"])
  
  #prep for building adjacency matrix
  d_sub<-d[!is.na(d$YIELD),]
  names<-unique(d_sub$GEOID)
  
  #clip county file to crop producing areas only and setup indexing by CNTY
  cnty<-cnty[cnty$GEOID %in% names,] %>% arrange(., GEOID) %>% mutate (.,CNTY=seq(1,nrow(.),1)) #order dataset by county, build group index that corresponds to adjacency matrix
  
  ##Build the relational matrix for areas of interest
  neighbors<-as(cnty,"Spatial") %>% poly2nb(., queen=F)#convert county sf to spatial polygons for inla relation matrix & create neighbors list from polygon object (neighbors share one or points at boundary)
  #temp <- poly2nb(shp, queen=T)#neighbors must share more than one point at boundary
  H.adj <- nb2mat(neighbors, style ="B", zero.policy=TRUE ) #convert to a sparse matrix to reduce memory (neighbor list to binary coded neighbor weights matrix)
  H.adj <-as(H.adj, "dgTMatrix") #convert to a sparse style matrix
  saveRDS(H.adj, paste0(wd,"/Data/H.adj.",crop,".rds"))
  
 
  ###add indexes and set up vars for log-linear model with rw (by binning) and quadratic (by squaring) terms
  d_sub<-left_join(d_sub,st_set_geometry(cnty[,c("GEOID","CNTY")], NULL),by="GEOID")#join new index to crop dataset
  d_sub<-d_sub %>% arrange(.,CNTY,YEAR) #Ordering dataset by county then time f
  
  d_sub$AERCODE.id<-d_sub$AERCODE
  d_sub$YEAR.id<-d_sub$Yr  #so index starts at 1
  d_sub$YEAR.id2<-(d_sub$Yr )^2
  d_sub$AERCODE.id2<-d_sub$AERCODE
  
  #Scale data and bin for rw1 evaluation
   d_sub<- d_sub %>% mutate_at(ivars, ~scale(.)) %>%
    mutate_at(c("TP","GDD","SDD", "SOIL", "PERC_IRR", "ACRES"), ~scale(.)) %>%
    mutate_at(ivars, ~inla.group(., n = 20, method = "quantile")) %>%   
    mutate_at(c("TP","GDD","SDD", "SOIL"), ~inla.group(., n = 20, method = "quantile")) 
  
  d_sub$PERC_IRR[is.na(d_sub$PERC_IRR)]<-0
  
  #Yield transformed to log(Yield)
  d_sub$YIELD<-log(d_sub$YIELD)
  
  
  saveRDS(d_sub,paste0(wd,"/Data/",savename,".rds"))
  
}


#################################################
###function to prep trimmed data for modeling (used for sensitivity testing only)
################################################
prep.data.trim<-function(file, savename, crop, projection){
  
  dat<-readRDS(paste0(wd,"/Data/",file,".RDS")) #read in the data
  dat<-dat %>% mutate(., Yr = as.numeric(as.character(YEAR))) %>% mutate(., Yr= Yr - min(Yr)) 
  dat$GEOID<-as.character(dat$GEOID)
  
  cnty<-readRDS(paste0(wd,"/Data/county.RDS"))
  cnty<-st_as_sf(cnty) %>% dplyr::select(.,c(5,11)) %>% st_transform(., projection)
  
  cnty_AER<-readRDS(paste0(wd,"/Data/cnty_AER_SA.RDS"))%>% st_transform(., projection)

  
  #join to add county and AER ids to dataset
  d<-left_join(cnty,dat, by="GEOID")
  d<-left_join(d,st_set_geometry(cnty_AER[,c(1,2)], NULL), by='GEOID')
  
  
  #prep for building adjacency matrix, dropping counties for which there are not at least 3 years of data and Percent Acreage in crop of interest is at least 5%
  d_sub<-d[!is.na(d$YIELD),]
  d_sub<-d_sub %>% group_by(GEOID) %>% mutate(.,keep=(n())) %>% ungroup()
  d_sub<- d_sub %>% filter(., keep>2) %>% dplyr::select(., -keep)
  d_sub<-d_sub %>% filter(., ACRES>=5) 
  names<-unique(d_sub$GEOID)
  
  #clip county file to crop producing areas only and setup indexing by CNTY
  cnty<-cnty[cnty$GEOID %in% names,] %>% arrange(., GEOID) %>% mutate (.,CNTY=seq(1,nrow(.),1)) #order dataset by county, build group index that corresponds to adjacency matrix
  
  ##Build the relational matrix for areas of interest
  neighbors<-as(cnty,"Spatial") %>% poly2nb(., queen=F)#convert county sf to spatial polygons for inla relation matrix & create neighbors list from polygon object (neighbors share one or points at boundary)
  #temp <- poly2nb(shp, queen=T)#neighbors must share more than one point at boundary
  H.adj <- nb2mat(neighbors, style ="B", zero.policy=TRUE ) #convert to a sparse matrix to reduce memory (neighbor list to binary coded neighbor weights matrix)
  H.adj <-as(H.adj, "dgTMatrix") #convert to a sparse style matrix
  saveRDS(H.adj, paste0(wd,"/Data/H.adj.trim",crop,".rds"))
  
 
  ###add indexes and set up vars for log-linear model with rw (by binning) and quadratic (by squaring) terms
  d_sub<-left_join(d_sub,st_set_geometry(cnty[,c("GEOID","CNTY")], NULL),by="GEOID")#join new index to crop dataset
  d_sub<-d_sub %>% arrange(.,CNTY,YEAR) #Ordering dataset by county then time for modeling as is done in the Ohio example in chapter 7 of r-inla book
  
  d_sub$AERCODE.id<-d_sub$AERCODE
  d_sub$YEAR.id<-d_sub$Yr  
  d_sub$YEAR.id2<-(d_sub$Yr )^2
  d_sub$AERCODE.id2<-d_sub$AERCODE
  
  
  d_sub<- d_sub %>% mutate_at(ivars, ~inla.group(., n = 20, method = "quantile")) %>%   #
    mutate_at(c("TP","GDD","SDD", "SOIL"), ~inla.group(., n = 20, method = "quantile")) #trimmed 7
  
  d_sub$PERC_IRR[is.na(d_sub$PERC_IRR)]<-0
  
  #Yield transformed to log(Yield)
  d_sub$YIELD<-log(d_sub$YIELD)
  
  
  saveRDS(d_sub,paste0(wd,"/Data/",savename,"trim.rds"))
  
}


####################################################################
###function to prep Data with added Landscape Indicator Quantiles 
#################################################################

data.quants<- function(file, savename, crop, projection){
  
  dat<-readRDS(paste0(wd,"/Data/",file,".RDS")) #read in the data
  dat<-dat %>% mutate(., Yr = as.numeric(as.character(YEAR))) %>% mutate(., Yr= Yr - min(Yr)) 
  dat$GEOID<-as.character(dat$GEOID)
  
  cnty<-readRDS(paste0(wd,"/Data/county.RDS"))
  cnty<-st_as_sf(cnty) %>% dplyr::select(.,c(5,11)) %>% st_transform(., projection)
  
  cnty_AER<-readRDS(paste0(wd,"/Data/cnty_AER_SA.RDS"))%>% st_transform(., projection)
  
  
  #join to add county and AER ids to dataset
  d<-left_join(cnty,dat, by="GEOID")
  d<-left_join(d,st_set_geometry(cnty_AER[,c(1,2)], NULL), by='GEOID')
  
  
  #prep for building adjacency matrix
  d_sub<-d[!is.na(d$YIELD),]
  names<-unique(d_sub$GEOID)
  
  #clip county file to crop producing areas only and setup indexing by CNTY
  cnty<-cnty[cnty$GEOID %in% names,] %>% arrange(., GEOID) %>% mutate (.,CNTY=seq(1,nrow(.),1)) #order dataset by county, build group index that corresponds to adjacency matrix
  
  #Build the relational matrix for areas of interest
  neighbors<-as(cnty,"Spatial") %>% poly2nb(., queen=F)#convert county sf to spatial polygons for inla relation matrix & create neighbors list from polygon object (neighbors share one or points at boundary)
  H.adj <- nb2mat(neighbors, style ="B", zero.policy=TRUE ) #convert to a sparse matrix to reduce memory (neighbor list to binary coded neighbor weights matrix)
  H.adj <-as(H.adj, "dgTMatrix") #convert to a sparse style matrix
  saveRDS(H.adj, paste0(wd,"/Data/H.adj.",crop,".rds"))
  
 
  #add indexes and set up vars for log-linear model with rw (by binning) and quadratic (by squaring) terms
  d_sub<-left_join(d_sub,st_set_geometry(cnty[,c("GEOID","CNTY")], NULL),by="GEOID")#join new index to crop dataset
  d_sub<-d_sub %>% arrange(.,CNTY,YEAR) #Ordering dataset by county then time for modeling as is done in the Ohio example in chapter 7 of r-inla book
  
  d_sub$AERCODE.id<-d_sub$AERCODE
  d_sub$YEAR.id<-d_sub$Yr  
  d_sub$YEAR.id2<-(d_sub$Yr )^2
  d_sub$AERCODE.id2<-d_sub$AERCODE
  
  #create quintile groupings for landscape metrics
  d_new<-d_sub %>% mutate_at(ivars, ~inla.group(., n = 5, method = "quantile")) %>%
    mutate_at(ivars, ~as.factor(.)) %>% dplyr::select(ivars) %>% 
    dplyr::rename_all(., ~tolower(gsub("LSM_","",.x))) %>% st_set_geometry(NULL)
  
  #combine original data with quintile data
  d_new<-dplyr::bind_cols(d_sub, d_new)
  
  #Scale data and bin for rw1 evaluation
  d_new<- d_new %>% mutate_at(ivars, ~scale(.)) %>%
    mutate_at(c("TP","GDD","SDD", "SOIL", "PERC_IRR", "ACRES"), ~scale(.)) %>%
    mutate_at(ivars, ~inla.group(., n = 20, method = "quantile")) %>%   
    mutate_at(c("TP","GDD","SDD", "SOIL"), ~inla.group(., n = 20, method = "quantile")) 
  
  d_new$PERC_IRR[is.na(d_new$PERC_IRR)]<-0
  
  #Yield transformed to log(Yield)
  d_new$YIELD<-log(d_new$YIELD)
  
  
  saveRDS(d_new,paste0(wd,"/Data/",savename,".rds"))
  
}
############################################
###function to build basic model formulas
###########################################
build.formula<-function(independent, crop, controls = "All", priors = c("pc", "default")){

  if (priors == "pc" & controls == "All"){       
    formula <- paste0('YIELD ~  1 + PERC_IRR + ACRES  + 
      f(SOIL, model="rw1", scale.model=T, hyper = list(theta = list(prior="pc.prec",param=c(v,0.01)))) +
      f(TP, model="rw1",scale.model=T, hyper = list (theta = list(prior="pc.prec",param=c(v,0.01)))) + 
      f(SDD, model="rw1",scale.model=T, hyper = list (theta = list(prior="pc.prec",param=c(v,0.01))))  + 
      f(GDD,model="rw1", scale.model=T, hyper = list (theta = list(prior="pc.prec",param=c(v,0.01)))) + 
      f(CNTY, model="bym2", graph=H.adj.',crop,', scale.model=TRUE, constr = TRUE,  
           hyper=list( phi =list(prior = "pc", param = c(phi.u, phi.alpha), inital=-3), 
                       prec =list(prior = "pc.prec", param = c(u,alpha), inital = 5))) +
      f(AERCODE.id, YEAR.id)+ f(AERCODE.id2, YEAR.id2)+ 
       f(',independent, ', model="rw1", scale.model=T, hyper = list (theta = list(prior="pc.prec",param=c(v,0.01))))')
    formula <- as.formula(noquote(gsub("[\n]","", formula)))
    }else{ 
     
  if (priors == "default" & controls == "All"){
    formula <- paste0('YIELD ~ 1 + PERC_IRR + ACRES + 
      f(SOIL, model="rw1", scale.model=T) +
      f(TP, model="rw1",scale.model=T) + 
      f(SDD,model="rw1",scale.model=T) + 
      f(GDD,model="rw1", scale.model=T) +
      f(CNTY, model="bym", graph=H.adj.', crop,', scale.model=TRUE) +
      f(AERCODE.id, YEAR.id)+ f(AERCODE.id2, YEAR.id2) +
      f(', independent, ', model="rw1", scale.model=T)')
    formula <- as.formula(noquote(gsub("[\n]","",formula)))
    }
    }
}

###############################################
###function to build formulas with constr (for sensitivity testing)
#############################################
build.formula.2<-function(independent, crop, controls = "All", priors = c("pc", "default")){
  
  if (priors == "pc" & controls == "All"){       
    formula <- paste0('YIELD ~  1 + PERC_IRR + ACRES  + 
      f(SOIL, model="rw1", scale.model=T, constr=TRUE, hyper = list(theta = list(prior="pc.prec",param=c(v,0.01)))) +
      f(TP, model="rw1",scale.model=T, constr=TRUE, hyper = list (theta = list(prior="pc.prec",param=c(v,0.01)))) + 
      f(SDD, model="rw1",scale.model=T, constr=TRUE,hyper = list (theta = list(prior="pc.prec",param=c(v,0.01))))  + 
      f(GDD,model="rw1", scale.model=T, constr=TRUE, hyper = list (theta = list(prior="pc.prec",param=c(v,0.01)))) + 
      f(CNTY, model="bym2", graph=H.adj.',crop,', scale.model=TRUE, constr = TRUE,  
           hyper=list( phi =list(prior = "pc", param = c(phi.u, phi.alpha), inital=-3), 
                       prec =list(prior = "pc.prec", param = c(u,alpha), inital = 5))) +
      f(AERCODE.id, YEAR.id)+ f(AERCODE.id2, YEAR.id2)+ 
       f(',independent, ', model="rw1", scale.model=T, constr=TRUE, hyper = list (theta = list(prior="pc.prec",param=c(v,0.01))))')
    formula <- as.formula(noquote(gsub("[\n]","", formula)))
  }else{ 
    
    if (priors == "default" & controls == "All"){
      formula <- paste0('YIELD ~ 1 + PERC_IRR + ACRES + 
      f(SOIL, model="rw1", scale.model=T) +
      f(TP, model="rw1",scale.model=T) + 
      f(SDD,model="rw1",scale.model=T) + 
      f(GDD,model="rw1", scale.model=T) +
      f(CNTY, model="bym", graph=H.adj.', crop,', scale.model=TRUE) +
      f(AERCODE.id, YEAR.id)+ f(AERCODE.id2, YEAR.id2) +
      f(', independent, ', model="rw1", scale.model=T)')
      formula <- as.formula(noquote(gsub("[\n]","",formula)))
    }
  }
}
 
##########################################
###function to build formulas with rw2 (sensitivity testing)
#########################################
build.formula.mod<-function(independent, crop, controls = "All", priors = c("pc", "default")){
  
  if (priors == "pc" & controls == "All"){       
    formula <- paste0('YIELD ~ -1 +  ACRES + PERC_IRR + 
      f(SOIL, model="rw2", scale.model=T, constr=TRUE, hyper = list(theta = list(prior="pc.prec",param=c(v,0.01)))) +
      f(TP, model="rw2",scale.model=T, constr=TRUE, hyper = list (theta = list(prior="pc.prec",param=c(v,0.01)))) + 
      f(SDD, model="rw2",scale.model=T, constr=TRUE, hyper = list (theta = list(prior="pc.prec",param=c(v,0.01))))  + 
      f(GDD,model="rw2", scale.model=T, constr=TRUE, hyper = list (theta = list(prior="pc.prec",param=c(v,0.01)))) + 
      f(CNTY, model="bym2", graph=H.adj.',crop,', scale.model=TRUE, constr = TRUE,  
           hyper=list( phi =list(prior = "pc", param = c(phi.u, phi.alpha), inital=-3), 
                       prec =list(prior = "pc.prec", param = c(u,alpha), inital = 5))) +
      f(AERCODE.id, YEAR.id)+ f(AERCODE.id2, YEAR.id2) +
       f(',independent, ', model="rw2", scale.model=T, constr=TRUE, hyper = list (theta = list(prior="pc.prec",param=c(v,0.01))))')
    formula <- as.formula(noquote(gsub("[\n]","", formula)))
  }else{ 
    
    if (priors == "default" & controls == "All"){
      formula <- paste0('YIELD ~ 1 + PERC_IRR + ACRES + SOIL +
      f(TP, model="rw1",scale.model=T) + 
      f(SDD,model="rw1",scale.model=T) + 
      f(GDD,model="rw1", scale.model=T) +
      f(CNTY, model="bym", graph=H.adj.', crop,', scale.model=TRUE) +
      f(AERCODE.id, YEAR.id)+ f(AERCODE.id2, YEAR.id2) +
      f(', independent, ', model="rw1", scale.model=T)')
      formula <- as.formula(noquote(gsub("[\n]","",formula)))
    }
  }
}


###########################################################
###function to build formulas without landscape metric (sensitivity testing)
########################################################
build.formula.nolm<-function(independent, crop, controls = "All", priors = c("pc", "default")){
  
  if (priors == "pc" & controls == "All"){       
    formula <- paste0('YIELD ~  1 + PERC_IRR + ACRES  + 
      f(SOIL, model="rw1", scale.model=T, hyper = list(theta = list(prior="pc.prec",param=c(v,0.01)))) +
      f(TP, model="rw1",scale.model=T, hyper = list (theta = list(prior="pc.prec",param=c(v,0.01)))) + 
      f(SDD, model="rw1",scale.model=T, hyper = list (theta = list(prior="pc.prec",param=c(v,0.01))))  + 
      f(GDD,model="rw1", scale.model=T, hyper = list (theta = list(prior="pc.prec",param=c(v,0.01)))) + 
      f(CNTY, model="bym2", graph=H.adj.',crop,', scale.model=TRUE, constr = TRUE,  
           hyper=list( phi =list(prior = "pc", param = c(phi.u, phi.alpha), inital=-3), 
                       prec =list(prior = "pc.prec", param = c(u,alpha), inital = 5))) +
      f(AERCODE.id, YEAR.id)+ f(AERCODE.id2, YEAR.id2))')
    formula <- as.formula(noquote(gsub("[\n]","", formula)))
  }else{ 
    
    if (priors == "default" & controls == "All"){
      formula <- paste0('YIELD ~ 1 + PERC_IRR + ACRES + 
      f(SOIL, model="rw1", scale.model=T) +
      f(TP, model="rw1",scale.model=T) + 
      f(SDD,model="rw1",scale.model=T) + 
      f(GDD,model="rw1", scale.model=T) +
      f(CNTY, model="bym", graph=H.adj.', crop,', scale.model=TRUE) +
      f(AERCODE.id, YEAR.id)+ f(AERCODE.id2, YEAR.id2) ')
      formula <- as.formula(noquote(gsub("[\n]","",formula)))
    }
  }
}
#################################################################################
###function to build formulas without landscape metrics and without climate/soil (sensitivity testing)
################################################################################
build.formula.nolmcs<-function(independent, crop, controls = "All", priors = c("pc", "default")){
  
  if (priors == "pc" & controls == "All"){       
    formula <- paste0('YIELD ~  1 + PERC_IRR + ACRES  + 
      f(CNTY, model="bym2", graph=H.adj.',crop,', scale.model=TRUE, constr = TRUE,  
           hyper=list( phi =list(prior = "pc", param = c(phi.u, phi.alpha), inital=-3), 
                       prec =list(prior = "pc.prec", param = c(u,alpha), inital = 5))) +
      f(AERCODE.id, YEAR.id)+ f(AERCODE.id2, YEAR.id2))')
    formula <- as.formula(noquote(gsub("[\n]","", formula)))
  }else{ 
    
    if (priors == "default" & controls == "All"){
      formula <- paste0('YIELD ~ 1 + PERC_IRR + ACRES + 
      f(CNTY, model="bym", graph=H.adj.', crop,', scale.model=TRUE) +
      f(AERCODE.id, YEAR.id)+ f(AERCODE.id2, YEAR.id2) ')
      formula <- as.formula(noquote(gsub("[\n]","",formula)))
    }
  }
}

###################################################
###function to build formulas without climate/soil (sensitivity testing)
###################################################
build.formula.nocs<-function(independent, crop, controls = "All", priors = c("pc", "default")){
  
  if (priors == "pc" & controls == "All"){       
    formula <- paste0('YIELD ~  1 + PERC_IRR + ACRES  + 
      f(CNTY, model="bym2", graph=H.adj.',crop,', scale.model=TRUE, constr = TRUE,  
           hyper=list( phi =list(prior = "pc", param = c(phi.u, phi.alpha), inital=-3), 
                       prec =list(prior = "pc.prec", param = c(u,alpha), inital = 5))) +
      f(AERCODE.id, YEAR.id)+ f(AERCODE.id2, YEAR.id2)) +
       f(',independent, ', model="rw2", scale.model=T, constr=TRUE, hyper = list (theta = list(prior="pc.prec",param=c(v,0.01))))')
    formula <- as.formula(noquote(gsub("[\n]","", formula)))
  }else{ 
    
    if (priors == "default" & controls == "All"){
      formula <- paste0('YIELD ~ 1 + PERC_IRR + ACRES + 
      f(CNTY, model="bym", graph=H.adj.', crop,', scale.model=TRUE) +
      f(AERCODE.id, YEAR.id)+ f(AERCODE.id2, YEAR.id2) +
      f(', independent, ', model="rw1", scale.model=T)')
      formula <- as.formula(noquote(gsub("[\n]","",formula)))
    }
  }
}


#############################################################################
###functions to build formulas with landscape interaction by quantile factor
############################################################################

#https://becarioprecario.bitbucket.io/inla-gitbook/ch-temporal.html
#create a seperable model using group function
#within-group effect (e.g. of SHDI in a single ED quantile group) is rw1
#between-group effect is iid

build.formula.facint1<-function(independent, crop, controls = "All", priors = c("pc", "default")){
  
  if (priors == "pc" & controls == "All"){       
    formula <- paste0('YIELD ~  1 + PERC_IRR + ACRES  + 
      f(SOIL, model="rw1", scale.model=T, hyper = list(theta = list(prior="pc.prec",param=c(v,0.01)))) +
      f(TP, model="rw1",scale.model=T, hyper = list (theta = list(prior="pc.prec",param=c(v,0.01)))) + 
      f(SDD, model="rw1",scale.model=T, hyper = list (theta = list(prior="pc.prec",param=c(v,0.01))))  + 
      f(GDD,model="rw1", scale.model=T, hyper = list (theta = list(prior="pc.prec",param=c(v,0.01)))) + 
      f(CNTY, model="bym2", graph=H.adj.',crop,', scale.model=TRUE, constr = TRUE,  
           hyper=list( phi =list(prior = "pc", param = c(phi.u, phi.alpha), inital=-3), 
                       prec =list(prior = "pc.prec", param = c(u,alpha), inital = 5))) +
      f(AERCODE.id, YEAR.id)+ f(AERCODE.id2, YEAR.id2)+ 
      f(',independent, ',  model="rw1", scale.model=T, group = as.numeric(ed_all), control.group = list(model ="iid"), hyper = list (theta = list(prior="pc.prec",param=c(v,0.01))))')
    formula <- as.formula(noquote(gsub("[\n]","", formula)))
  }else{ 
    
    if (priors == "default" & controls == "All"){
      formula <- paste0('YIELD ~ 1 + PERC_IRR + ACRES + 
      f(SOIL, model="rw1", scale.model=T) +
      f(TP, model="rw1",scale.model=T) + 
      f(SDD,model="rw1",scale.model=T) + 
      f(GDD,model="rw1", scale.model=T) +
      f(CNTY, model="bym", graph=H.adj.', crop,', scale.model=TRUE) +
      f(AERCODE.id, YEAR.id)+ f(AERCODE.id2, YEAR.id2) +
      f(', independent, ',  model="rw1", scale.model=T, group = as.numeric(ed_all), control.group = list(model ="iid"))')
      formula <- as.formula(noquote(gsub("[\n]","",formula)))
    }
  }
}

build.formula.facint2<-function(independent, crop, controls = "All", priors = c("pc", "default")){
  
  if (priors == "pc" & controls == "All"){       
    formula <- paste0('YIELD ~  1 + PERC_IRR + ACRES  + 
      f(SOIL, model="rw1", scale.model=T, hyper = list(theta = list(prior="pc.prec",param=c(v,0.01)))) +
      f(TP, model="rw1",scale.model=T, hyper = list (theta = list(prior="pc.prec",param=c(v,0.01)))) + 
      f(SDD, model="rw1",scale.model=T, hyper = list (theta = list(prior="pc.prec",param=c(v,0.01))))  + 
      f(GDD,model="rw1", scale.model=T, hyper = list (theta = list(prior="pc.prec",param=c(v,0.01)))) + 
      f(CNTY, model="bym2", graph=H.adj.',crop,', scale.model=TRUE, constr = TRUE,  
           hyper=list( phi =list(prior = "pc", param = c(phi.u, phi.alpha), inital=-3), 
                       prec =list(prior = "pc.prec", param = c(u,alpha), inital = 5))) +
      f(AERCODE.id, YEAR.id)+ f(AERCODE.id2, YEAR.id2)+ 
      f(',independent, ',  model="rw1", scale.model=T, group = as.numeric(area_mn_all), control.group = list(model ="iid"), hyper = list (theta = list(prior="pc.prec",param=c(v,0.01))))')
    formula <- as.formula(noquote(gsub("[\n]","", formula)))
  }else{ 
    
    if (priors == "default" & controls == "All"){
      formula <- paste0('YIELD ~ 1 + PERC_IRR + ACRES + 
      f(SOIL, model="rw1", scale.model=T) +
      f(TP, model="rw1",scale.model=T) + 
      f(SDD,model="rw1",scale.model=T) + 
      f(GDD,model="rw1", scale.model=T) +
      f(CNTY, model="bym", graph=H.adj.', crop,', scale.model=TRUE) +
      f(AERCODE.id, YEAR.id)+ f(AERCODE.id2, YEAR.id2) +
      f(', independent, ',  model="rw1", scale.model=T, group = as.numeric(area_mn_all), control.group = list(model ="iid"))')
      formula <- as.formula(noquote(gsub("[\n]","",formula)))
    }
  }
}
######################################## 
###function to run basic gaussian model
########################################
run.model<-function(formula, data, crop, priors ="pc", name){
  
  df<-readRDS(paste0(wd,"/Data/",data,".RDS")) #read in the data #define the dataset for the crop of interest
  h<-readRDS(paste0(wd,"/Data/H.adj.", crop, ".RDS"))
  assign(paste0("H.adj.", crop), h)
   v<-NA
   n<-NA
   Q<-NA
   u<-NA
   alpha<-NA
   phi.u<-NA
   phi.alpha<-NA
   phi.prior<-NA
  if (priors == "pc"){ #setup pc priors
    v=sd(df$YIELD,na.rm=TRUE)/0.31
    n=dim(df)[1] 
    Q = INLA:::inla.pc.bym.Q(h)
    Q = INLA:::inla.scale.model(Q, constr=list(A=matrix(1, 1, n), e=0))             
    u = 0.2/0.31
    alpha = 0.01
    phi.u = 0.5
    phi.alpha = 2/3 ## prob(phi < phi.u) = phi.alpha   
    phi.prior = INLA:::inla.pc.bym.phi(Q=Q, u= phi.u, alpha = phi.alpha)
  }
  OUT<-inla(formula, data=df, family = "gaussian",     #run the inla model
            control.predictor = list( compute=TRUE), 
            control.compute=list(dic=TRUE, cpo=TRUE), 
            control.fixed = list(prec = 0.0000001)) 
  saveRDS(OUT, file=paste0(wd,"/Output/",name ,".rds")) #save the model output
 
} 

################################################################################
###function to return a summary table of model results and if desired diagnotics
#################################################################################
model.results<-function(model_results_name, data="d", crop="corn", diagnostics = TRUE){
  results<-readRDS(paste0(wd,"/Output/",model_results_name,".rds"))
  data<-readRDS(paste0(wd,"/Data/",data,".RDS"))
 # st_geometry(data)<-NULL
  summary<-summary(results)
  if (diagnostics == TRUE){ #if diagnostics turned on return the dic, cpo, and pit
    #R_INLA Diagnostics
    DIC<-(results$dic$dic)
    CPO<-(sum(log(results$cpo$cpo), na.rm=T)) 
    h<-ggplot(as.data.frame(results$cpo$pit), aes(x=results$cpo$pit)) + geom_histogram() + ggtitle( "Model PIT")
    #PPC Distribution Checks
    p_val<-c()
    n<-length(data$YIELD)
    for(i in (1:n)){
      p_val[i]<-inla.pmarginal(q=data$YIELD[i],
                               marginal=results$marginals.fitted.values[[i]])
    }
    data<-cbind(data,results$summary.fitted.values$mean)
    colnames(data)[(dim(data)[2]-1)]<-c("FittedVals")
    p1<-ggplot(data[!is.na(data$YIELD),], aes(x=YIELD,y=FittedVals)) + geom_point(size=0.2, color="blue") +
      ggtitle("Scatter Plot of Predicted and Observed Values")  + xlab("Observed Value") + ylab("Mean of Post. Pred. Distr.")
    
    p2<-ggplot(as.data.frame(p_val), aes(x=p_val)) + geom_histogram() + ggtitle("Posterior Predictive p-values")
    
    #PPC Summary Metrics
    sq_dif<-(data$YIELD-data$FittedVals)^2
    MSE<-1/n*(sum(sq_dif,na.rm=T))
    
    pred_res2<-(data$FittedVals[!is.na(data$YIELD)] - mean(data$YIELD, na.rm=T)) ^2
    obs_res2<-(data$YIELD[!is.na(data$YIELD)] - mean(data$YIELD, na.rm=T))^2
    R2<-sum(pred_res2, na.rm=T)/sum(obs_res2, na.rm=T)
    
    #OUPUTS
    summtable<-as.data.frame(summary$fixed[,1:5])
    summtable.pt2<-as.data.frame(summary$hyperpar[,1:5])
    summtable<-rbind(summtable,summtable.pt2)
    summtable<-kable(summtable, caption="Summary Table of Model Estimates")
    
    Diagnostics<-as.data.frame(t(c(DIC,CPO,MSE,R2)))
    colnames(Diagnostics)<-c("DIC","CPO","MSE","R2")
    Diagnostics<-kable(Diagnostics, caption="Model Diagnostic Metrics")
    
    figure<-ggarrange(h,p1,p2,ncol=3,nrow=1)
    annotate_figure(figure, top= text_grob("Model Diagnostic Plots"))
    
    my_list<-list(summtable,  figure, Diagnostics)
    return(my_list)
  }  
  my_list<-list(summary)
  return(my_list)
}


##############################################
###function to return maps of spatial effects
##############################################
plot.spatialeffect<-function(model_results_name, data, crop, type ="BYM", scale="CNTY"){
  data<-readRDS(paste0(wd,"/Data/",data,".RDS"))
  results<-readRDS(paste0(wd,"/Output/",model_results_name,".rds"))
  
  if (type =="BYM"){ #Spatial effect for bym models
    if (scale =="CNTY"){ #County spatial effects
      Nareas<-length(unique(data$CNTY)) #number of unique spatial locations in dataset
      asr<-results$summary.random$CNTY[1:Nareas,c(1,2,3)]#extract the area-specific residuals
      
      #Create spatial dataframe with information for map and plot
      map.asr<-left_join(data[,c("CNTY","GEOID","geometry")],asr,by=c("CNTY"="ID"))
      colnames(map.asr)<-c("CNTY","GEOID","Mean","Standard Deviation","geometry")
      
      saveRDS(map.asr,paste0(wd,"/output/",model_results_name,"_asr.rds"))
      p1<-plot(map.asr[,3], main= paste("Mean of Area Specific Residuals for"), lwd=0.5, key.pos=1)
      p2<-plot(map.asr[,4], main=paste(model_results_name,"Standard Deviation of Area Specific Residuals", crop), lwd=0.5, key.pos=1)
      
      
      my_list<-list(p1,p2)
      
      return (my_list)
    }
    
  }
  
}


#######################################################
###function to plot smooth effects of independent vars
##########################################################
plot.smootheffect<-function(model_results_name, var="SDI"){
  
  results<-readRDS(paste0(wd,"/Output/",model_results_name,".rds"))
  
    p2<-ggplot(data=results$summary.random[[var]][ ,c(1,4,5,6)], aes (x=ID, y=`0.5quant`)) + geom_point() +
      geom_ribbon(aes(x=ID,ymin = `0.025quant`, ymax = `0.975quant`, fill="red"), alpha=0.25) + 
      theme(legend.position="none")+ xlab(var) + ylab("Estimated Effect on log(Yield)")
    return(p2)
    
}


#########################################################################################################
###function to create summary plots combining all three crops effect estimates for each landscape metric
#########################################################################################################
plot.lm.summ<-function(metrics, crops, y.min, y.max, x.min, x.max){
  
  
  metric<-toupper(metrics)
  metrics.labels <-c("Mean Patch Area","Contagion","Edge Density","Largest Patch Index",
                     "Richness","Shannon Diversity Index","Shannon Evenness Index","Simpson Diversity Index","Simpson Evenness Index","Percent Natural Cover")
  crops.labels<-crops
  crops.labels[crops.labels == "wwheat"] <- "Winter wheat"
  crops.labels[crops.labels == "corn"] <- "Corn"
  crops.labels[crops.labels == "soy"] <- "Soy"
  
  
  my.plots<-c()
  
  foreach (i=1:length(metrics))%do% {
    
    obs<-NA
    
    foreach (j=1:length(crops))%do% {
      #pull the model output by metric for ag diversity
      if (file.exists(paste0(wd,"/Output/", crops[j], "_",metrics[i],".rds"))){
        results_1<-readRDS(paste0(wd,"/Output/", crops[j], "_",metrics[i],".rds"))
        results_1<-get(metric[i],results_1$summary.random)[ ,c(1,4,5,6)] #pull the results for the effect of the metric for the crop
        results_1$Crop<-crops.labels[j]#add crop column
        obs<-c(obs, list(results_1))
      }
    }
    
    
    #combine data for a single plot
    library(data.table)
    
    results_tot<-as.data.frame(rbindlist(obs[2:4]))
    mean_lm<-results_tot %>% group_by (Crop) %>% summarise_at('ID', mean)
    sd_lm<-results_tot %>% group_by (Crop) %>% summarise_at('ID', sd)
    sd_lm$ID<-sd_lm$ID+mean_lm$ID
    sdlow_lm<-results_tot %>% group_by (Crop) %>% summarise_at('ID', sd)
    sdlow_lm$ID<-mean_lm$ID-sdlow_lm$ID
    
   
    #build a plot
    library("viridis")
    
    p1<-ggplot(data=results_tot) + geom_line(aes (x=ID, y=`0.5quant`, group = Crop, col=Crop), size=0.5) +
      xlim(x.min[i], x.max[i]) + ylim(y.min,y.max) +
      xlab(paste0(metrics.labels[i])) + ylab("Effect on log (Yield)") +   
      scale_x_continuous(breaks=c(-3,-2,-1, 0, 1,2,3), limits=c(x.min[i], x.max[i])) + #this overwrites the min-max limit which we want to be the same for all plots
      geom_ribbon(aes(x=ID, ymin = `0.025quant`, ymax = `0.975quant`, group=Crop, fill=Crop), alpha=0.25) +
      scale_color_viridis(discrete = TRUE, option = "D")+
      scale_fill_viridis(discrete = TRUE, option = "D") +
      # geom_vline(aes(xintercept=ID, col=Crop), mean_lm) +
      # geom_vline(aes(xintercept=ID, col=Crop), sd_lm, linetype="dashed") +
      # geom_vline(aes(xintercept=ID, col=Crop), sdlow_lm, linetype="dashed") +
      theme_light() + theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        # panel.grid.major = element_line(colour = "grey50"),
        #   panel.grid.minor = element_line(colour = "grey50"),
        # panel.background = element_rect(fill = NA),
        panel.ontop = FALSE
        )    
    
    my.plots[[i]]<-p1
    
  }
  return(my.plots)
}


######################################################################################################################
###function to create simple summary plots combining all climate response curves for each crop-landscape metric model
#######################################################################################################################
plot.cl.summ.simple<-function(metrics, crops, y.min, y.max, x.min, x.max){
  
  metrics.labels<-toupper(metrics)
  crops.labels<-crops
  crops.labels[crops.labels == "wwheat"] <- "Winter wheat"
  crops.labels[crops.labels == "corn"] <- "Corn"
  crops.labels[crops.labels == "soy"] <- "Soy"
  climate.labels<-c("TP","GDD","SDD","SOIL")
  climate.labels.full<-c("Total Precipitation","Growing Degree Days","Stress Degree Days","Soil Suitability")
  
  
  my.plots<-c()
  
  foreach (i=1:length(crops))%do% {
    
    foreach(j=1:length(climate.labels))%do%{
      
      obs<-NA
      
      foreach (k=1:length(metrics))%do% {
        #pull the model output by metric for ag diversity
        if (file.exists(paste0(wd,"/Output/", crops[i], "_",metrics[k],".rds"))){
          results_1<-readRDS(paste0(wd,"/Output/", crops[i], "_",metrics[k],".rds"))
          results_1<-get(climate.labels[j],results_1$summary.random)[ ,c(1,4,5,6)] #pull the results for the effect of the metric for the crop
          results_1$Metric<-metrics.labels[k]#add metric column
          results_1$Climate<-climate.labels[j]
          results_1$Crop<-crops.labels[i]
          obs<-c(obs, list(results_1))
        }
      }
      
      
      #combine data for a single plot
      library(data.table)
      
      results_tot<-as.data.frame(rbindlist(obs[2:(length(metrics)+1)])) #combine all landscape metric models
      
      if (j==1){
        results_tot_cl <- results_tot
      }else{
      results_tot_cl<-as.data.frame(rbind(results_tot_cl, results_tot))
      }
    }
    
    if (i==1){
      results_tot_cl_crop <-results_tot_cl    
    }else{
        results_tot_cl_crop<-as.data.frame(rbind(results_tot_cl_crop,results_tot_cl))
    }
  }
  
  #build plots
  foreach(j=1:length(climate.labels))%do%{
      
      p1<-results_tot_cl_crop %>% filter(Climate==climate.labels[j]) %>% 
        ggplot() + geom_line(aes (x=ID, y=`0.5quant`, group = interaction(Metric, Crop), col= Crop), size=0.5) +
        xlim(x.min[j], x.max[j]) + ylim(y.min[j],y.max[j]) +
        xlab(paste0(climate.labels.full[j])) + ylab("Effect on log (Yield)") +  
        geom_ribbon(aes(x=ID, ymin = `0.025quant`, ymax = `0.975quant`, group=interaction(Metric, Crop), fill=Crop), alpha=0.025) +
        scale_color_viridis(discrete = TRUE, option = "D")+
        scale_fill_viridis(discrete = TRUE, option = "D") +
        theme_light() + 
        theme(panel.grid.major.x = element_blank(),
               panel.grid.minor.x = element_blank(),
                panel.ontop = FALSE)
      
      my.plots[[j]]<-p1
    }
   
  return(my.plots)
}

###############################################################################################################
###function to create summary plots combining all climate response curves for each crop-landscape metric model
###############################################################################################################
plot.lm.summ.cl<-function(metrics, crops, y.min, y.max, x.min, x.max){
  
  metrics.labels<-toupper(metrics)
  crops.labels<-crops
  crops.labels[crops.labels == "wwheat"] <- "Winter wheat"
  crops.labels[crops.labels == "corn"] <- "Corn"
  crops.labels[crops.labels == "soy"] <- "Soy"
  climate.labels<-c("TP","GDD","SDD","SOIL")
  climate.labels.full<-c("Total Precipitation","Growing Degree Days","Stress Degree Days","Soil Suitability")
  
  all.plots<-c()
  my.plots<-c()
  
  foreach (i=1:length(crops))%do% {
    
    foreach(j=1:length(climate.labels))%do%{
    
        obs<-NA
    
        foreach (k=1:length(metrics))%do% {
        #pull the model output by metric for ag diversity
            if (file.exists(paste0(wd,"/Output/", crops[i], "_",metrics[k],".rds"))){
              results_1<-readRDS(paste0(wd,"/Output/", crops[i], "_",metrics[k],".rds"))
              results_1<-get(climate.labels[j],results_1$summary.random)[ ,c(1,4,5,6)] #pull the results for the effect of the metric for the crop
              results_1$Metric<-metrics.labels[k]#add metric column
              results_1$Climate<-climate.labels[j]
              results_1$Crop<-crops.labels[i]
              obs<-c(obs, list(results_1))
            }
        }
    
    
    #combine data for a single plot
    library(data.table)
    
    results_tot<-as.data.frame(rbindlist(obs[2:(length(metrics)+1)]))
    
    #build a plot
    p1<-ggplot(data=results_tot) + geom_line(aes (x=ID, y=`0.5quant`, group = Metric, col= Metric), size=0.5) +
      xlim(x.min[j,i], x.max[j,i]) + ylim(y.min[j,i],y.max[j,i]) +
      xlab(paste0(climate.labels.full[j])) + ylab("Effect on log (Yield)") +  ggtitle (paste0(crops.labels[i])) +
      geom_ribbon(aes(x=ID, ymin = `0.025quant`, ymax = `0.975quant`, group=Metric), alpha=0.025) +
      theme_classic()
    
    my.plots[[j]]<-p1
    
    }
    all.plots[[i]]<-my.plots
    }
  return(all.plots)
}

################################################################################
###function to return a summary table of model fit statistics for models
#############################################################################
model.results.summ<-function(data, crop="corn", metrics=metrics){
  
  library(flextable)
  Diagnostics<-data.frame(matrix(ncol=6, nrow=length(metrics)))
  colnames(Diagnostics)<-c("Metric", "DIC","CPO","MSE","Mean P Value", "R2")
  
  foreach (i=1:length(metrics))%do% {
    
    results<-readRDS(paste0(wd,"/Output/",crop, "_", metrics[i],".rds"))
    d<-readRDS(paste0(wd,"/Data/",data,".RDS"))
    
    # return the dic, cpo
    #R_INLA Diagnostics
    Diagnostics$DIC[i]<-(results$dic$dic)
    Diagnostics$CPO[i]<-(sum(log(results$cpo$cpo), na.rm=T)) 
    
    #PPC Distribution Checks
    p_val<-c()
    n<-length(d$YIELD)
    for(j in (1:n)){
      p_val[j]<-inla.pmarginal(q=d$YIELD[j],
                               marginal=results$marginals.fitted.values[[j]])
    }
    Diagnostics$`Mean P Value`[i]<-mean(p_val, na.rm=T)
    
    d<-cbind(d,results$summary.fitted.values$mean)
    colnames(d)[(dim(d)[2]-1)]<-c("FittedVals")
    
    
    #PPC Summary Metrics
    sq_dif<-(d$YIELD-d$FittedVals)^2
    Diagnostics$MSE[i]<-1/n*(sum(sq_dif,na.rm=T))
    
    pred_res2<-(d$FittedVals[!is.na(d$YIELD)] - mean(d$YIELD, na.rm=T)) ^2
    obs_res2<-(d$YIELD[!is.na(d$YIELD)] - mean(d$YIELD, na.rm=T))^2
    Diagnostics$R2[i]<-sum(pred_res2, na.rm=T)/sum(obs_res2, na.rm=T)
    Diagnostics$Metric[i]<-metrics[i]
    #OUPUTS
    
  }
  
  
  Diag<-regulartable(Diagnostics)
  print(Diag, preview="docx")
  return(Diagnostics)
}


##########################################################
###function to create an interaction summary plot
###############################################################
plot.interaction<-function(metrics, crops, y.min, y.max, x.min, x.max){
  
  metrics.labels<-toupper(metrics)
  crops.labels<-crops
  crops.labels[crops.labels == "wwheat"] <- "Winter wheat"
  crops.labels[crops.labels == "corn"] <- "Corn"
  crops.labels[crops.labels == "soy"] <- "Soy"
  composition.labels<-c("LSM_SHDI_ALL","LSM_SHDI_ALL","LSM_RICH_ALL","LSM_RICH_ALL")
  configuration.labels<-c("Edge Density","Mean Area","Edge Density","Mean Area")
  
  
  my.plots<-c()
  
  foreach (i=1:length(crops))%do% {
    
    # foreach(j=1:length(climate.labels))%do%{
      
      obs<-NA
      
      foreach (k=1:length(metrics))%do% {
        #pull the model output by metric for ag diversity
        if (file.exists(paste0(wd,"/Output/", crops[i], "_",metrics[k],"quants.rds"))){
          results_1<-readRDS(paste0(wd,"/Output/", crops[i], "_",metrics[k],"quants.rds"))
          results_1<-get(composition.labels[k],results_1$summary.random)[ ,c(1,4,5,6)] #pull the results for the effect of the metric for the crop
          results_1$Metric<-metrics.labels[k]#add metric column
          results_1$Composition<-composition.labels[k]
          results_1$Configuration<-configuration.labels[k]
          results_1$Crop<-crops.labels[i]
          quintiles<-c("1st Quintile","2nd Quintile", "3rd Quintile","4th Quintile","5th Quintile")
            if (configuration.labels[k]=="Mean Area"){ #Since complexity increases with decreasing Mean Area reverse the order of quintile labels
              quintiles <- rev(quintiles)
            }
            if (composition.labels[k]=="LSM_RICH_ALL" & crops[i]=="soy"){
          results_1$Quintile<-rep(quintiles,each=19) 
            }else{
              results_1$Quintile<-rep(quintiles,each=20)
            }
          obs<-c(obs, list(results_1))
        }
      }
      
      
      #combine data for a single plot
      library(data.table)
      
      results_tot<-as.data.frame(rbindlist(obs[2:(length(metrics)+1)])) #combine all landscape metric models
      
    
    p1<-results_tot %>% filter(Crop==crops.labels[i]) %>% 
      ggplot() + geom_line(aes (x=ID, y=`0.5quant`, group = interaction(Quintile, Metric), col= Quintile), size=0.5) +
      xlim(x.min, x.max) + ylim(y.min,y.max) +
      xlab(paste0('Composition')) + ylab("Effect on log (Yield)") +  
      ggtitle(paste0(crops.labels[i])) +
      geom_ribbon(aes(x=ID, ymin = `0.025quant`, ymax = `0.975quant`, group=interaction(Quintile, Metric), fill=Quintile), alpha=0.05) +
      scale_color_viridis(discrete = TRUE, option = "D", name="Configuration")+
      scale_fill_viridis(discrete = TRUE, option = "D", name="Configuration") +
      theme_light() + theme(panel.grid.major.x = element_blank(),
                            panel.grid.minor.x = element_blank(),
                            panel.ontop = FALSE) 
    
   my.plots[[i]]<-p1 
  
  }
  return(my.plots)
}

