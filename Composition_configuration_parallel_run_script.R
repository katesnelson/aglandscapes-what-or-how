library(pacman)

p_load(tidyverse, INLA, sf, spdep, rgeos, knitr, ggpubr,doParallel, foreach,tidyselect, viridis)


wd<-getwd()

source(paste0(wd,'/Scripts/SA_functions.R'))

############################################
### Run models for each landscape metric ###
############################################

   ivars<-c("LSM_AREA_MN_ALL", "LSM_CONTAG_ALL","LSM_ED_ALL","LSM_LPI_ALL",
            "LSM_RICH_ALL","LSM_SHDI_ALL","LSM_SHEI_ALL","LSM_SIDI_ALL","LSM_SIEI_ALL", "PNC")
   datafiles<-c("corn_panel", "soy_panel","wwheat_panel")
   crops<-c("corn","soy","wwheat")
   
   foreach (j=1:length(datafiles))%do%{
   prep.data.std(file=datafiles[j], savename=paste0(crops[j],"std"), crop=crops[j], 
                 projection="+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m no_defs") #102003 epsg code for USA albers equal area conic
   }
   
   
   cl <- makeCluster(length(ivars))
   registerDoParallel(cl)
   foreach(i=1:length(ivars), .packages=c("sf", "dplyr","INLA","doParallel","foreach")) %dopar% {
   tryCatch({
       foreach (j=1:length(crops)) %do%{
         
         f1<-build.formula(ivars[i], crops[j], controls = "All", priors = c("default"))
         m1<- run.model(f1, data=paste0(crops[j],"std"), crop=crops[j], priors ="default", name=paste0(crops[j],"_",ivars[i])) 
         
      }
     }, error=function(cond) {
           message(paste("Encountered some issue for", ivars[i], "and", crops[j])) #report that there is an error and where it occured, not working, add to fcn?
           
         })
   }
   stopCluster(cl)


###################################
### Return model fit statistics ###
###################################
   
   model.results.summ(data="cornfinal", crop="corn", metrics)
   model.results.summ(data="soyfinal", crop="soy", metrics)
   model.results.summ(data="wwheatfinal", crop="wwheat", metrics)

###############################
### Sensitivity test models ###
###############################
   
   #Run model without diversity
   
         f1<-build.formula.nolm(ivars[1], crops[1], controls = "All", priors = c("default"))
         m1<- run.model(f1, data=paste0(crops[1],"final"), crop=crops[1], priors ="default", name=paste0(crops[1],"_",ivars[1],"nolm")) #set so that if pc prior model fails it will try to run default priors
         model.results.summ.trim(data="cornfinal", crop="corn",metrics[1])
         
   #Run model without climate/soil and landscape
         
         f1<-build.formula.nolmcs(ivars[1], crops[1], controls = "All", priors = c("default"))
         m1<- run.model(f1, data=paste0(crops[1],"final"), crop=crops[1], priors ="default", name=paste0(crops[1],"_",ivars[1],"nolmcs")) #set so that if pc prior model fails it will try to run default priors
         model.results.summ.trim(data="cornfinal", crop="corn",metrics[1])
   
   #Run model without climate/soil
         
         f1<-build.formula.nocs(ivars[1], crops[1], controls = "All", priors = c("default"))
         m1<- run.model(f1, data=paste0(crops[1],"final"), crop=crops[1], priors ="default", name=paste0(crops[1],"_",ivars[1],"nocs")) #set so that if pc prior model fails it will try to run default priors
         model.results.summ.trim(data="cornfinal", crop="corn",metrics[1])
      
 
#################################
### Extract summary measures ### 
################################
   
#Get max impact (median effect), total change (conservative, using confidence bounds), and -1 to +1 sd change (median effect) summary measures      
      
      #RICH
      p1[[5]]$data %>% filter(Crop =="Winter wheat") %>% dplyr::select(starts_with("0.5")) %>% max()
      p1[[5]]$data %>% filter(Crop =="Winter wheat") %>% dplyr::select(starts_with("0.025")| starts_with("0.975")) %>% 
        mutate_at(., .vars=c(1), ~max(.)) %>% mutate_at(., .vars=c(2), ~min(.)) %>% mutate(., change=.[,1]-.[,2]) %>% .[1, "change"]
      p1[[5]]$data %>% filter(Crop =="Winter wheat") %>% dplyr::select(starts_with("0.5")| ID) %>% 
        mutate(top = .[near(ID,1, tol=0.05),1]) %>% mutate(bottom = .[near(ID,-1, tol=0.05),1]) %>% mutate(., change=.[,3]-.[,4]) %>% .[1, "change"]
      
      p1[[5]]$data %>% filter(Crop =="Corn") %>% dplyr::select(starts_with("0.5")) %>% max()
      p1[[5]]$data %>% filter(Crop =="Corn") %>% dplyr::select(starts_with("0.025")| starts_with("0.975")) %>% 
        mutate_at(., .vars=c(1), ~max(.)) %>% mutate_at(., .vars=c(2), ~min(.)) %>% mutate(., change=.[,1]-.[,2]) %>% .[1, "change"]
      p1[[5]]$data %>% filter(Crop =="Corn") %>% dplyr::select(starts_with("0.5")| ID) %>% 
        mutate(top = .[near(ID,1, tol=0.1),1]) %>% mutate(bottom = .[near(ID,-1, tol=0.05),1]) %>% mutate(., change=.[,3]-.[,4]) %>% .[1, "change"]
      
      p1[[5]]$data %>% filter(Crop =="Soy") %>% dplyr::select(starts_with("0.5")) %>% max()
      p1[[5]]$data %>% filter(Crop =="Soy") %>% dplyr::select(starts_with("0.025")| starts_with("0.975")) %>% 
        mutate_at(., .vars=c(1), ~max(.)) %>% mutate_at(., .vars=c(2), ~min(.)) %>% mutate(., change=.[,1]-.[,2]) %>% .[1, "change"]
      p1[[5]]$data %>% filter(Crop =="Soy") %>% dplyr::select(starts_with("0.5")| ID) %>% 
        mutate(top = .[near(ID,1, tol=0.2),1]) %>% mutate(bottom = .[near(ID,-1, tol=0.05),1]) %>% mutate(., change=.[,3]-.[,4]) %>% .[1, "change"]
      
      #SHDI
      p1[[6]]$data %>% filter(Crop =="Winter wheat") %>% dplyr::select(starts_with("0.5")) %>% max()
      p1[[6]]$data %>% filter(Crop =="Winter wheat") %>% dplyr::select(starts_with("0.025")| starts_with("0.975")) %>% 
        mutate_at(., .vars=c(1), ~max(.)) %>% mutate_at(., .vars=c(2), ~min(.)) %>% mutate(., change=.[,1]-.[,2]) %>% .[1, "change"]
      p1[[6]]$data %>% filter(Crop =="Winter wheat") %>% dplyr::select(starts_with("0.5")| ID) %>% 
        mutate(top = .[near(ID,1, tol=0.1),1]) %>% mutate(bottom = .[near(ID,-1, tol=0.05),1]) %>% mutate(., change=.[,3]-.[,4]) %>% .[1, "change"]
      
      p1[[6]]$data %>% filter(Crop =="Corn") %>% dplyr::select(starts_with("0.5")) %>% max()
      p1[[6]]$data %>% filter(Crop =="Corn") %>% dplyr::select(starts_with("0.025")| starts_with("0.975")) %>% 
        mutate_at(., .vars=c(1), ~max(.)) %>% mutate_at(., .vars=c(2), ~min(.)) %>% mutate(., change=.[,1]-.[,2]) %>% .[1, "change"]
      p1[[6]]$data %>% filter(Crop =="Corn") %>% dplyr::select(starts_with("0.5")| ID) %>% 
        mutate(top = .[near(ID,1, tol=0.1),1]) %>% mutate(bottom = .[near(ID,-1, tol=0.05),1]) %>% mutate(., change=.[,3]-.[,4]) %>% .[1, "change"]
      
      p1[[6]]$data %>% filter(Crop =="Soy") %>% dplyr::select(starts_with("0.5")) %>% max()
      p1[[6]]$data %>% filter(Crop =="Soy") %>% dplyr::select(starts_with("0.025")| starts_with("0.975")) %>% 
        mutate_at(., .vars=c(1), ~max(.)) %>% mutate_at(., .vars=c(2), ~min(.)) %>% mutate(., change=.[,1]-.[,2]) %>% .[1, "change"]
      p1[[6]]$data %>% filter(Crop =="Soy") %>% dplyr::select(starts_with("0.5")| ID) %>% 
        mutate(top = .[near(ID,1, tol=0.1),1]) %>% mutate(bottom = .[near(ID,-1, tol=0.05),1]) %>% mutate(., change=.[,3]-.[,4]) %>% .[1, "change"]
      
      #SHEI
      p1[[7]]$data %>% filter(Crop =="Winter wheat") %>% dplyr::select(starts_with("0.5")) %>% max()
      p1[[7]]$data %>% filter(Crop =="Winter wheat") %>% dplyr::select(starts_with("0.025")| starts_with("0.975")) %>% 
        mutate_at(., .vars=c(1), ~max(.)) %>% mutate_at(., .vars=c(2), ~min(.)) %>% mutate(., change=.[,1]-.[,2]) %>% .[1, "change"]
      p1[[7]]$data %>% filter(Crop =="Winter wheat") %>% dplyr::select(starts_with("0.5")| ID) %>% 
        mutate(top = .[near(ID,1, tol=0.1),1]) %>% mutate(bottom = .[near(ID,-1, tol=0.05),1]) %>% mutate(., change=.[,3]-.[,4]) %>% .[1, "change"]
      
      p1[[7]]$data %>% filter(Crop =="Corn") %>% dplyr::select(starts_with("0.5")) %>% max()
      p1[[7]]$data %>% filter(Crop =="Corn") %>% dplyr::select(starts_with("0.025")| starts_with("0.975")) %>% 
        mutate_at(., .vars=c(1), ~max(.)) %>% mutate_at(., .vars=c(2), ~min(.)) %>% mutate(., change=.[,1]-.[,2]) %>% .[1, "change"]
      p1[[7]]$data %>% filter(Crop =="Corn") %>% dplyr::select(starts_with("0.5")| ID) %>% 
        mutate(top = .[near(ID,1, tol=0.1),1]) %>% mutate(bottom = .[near(ID,-1, tol=0.05),1]) %>% mutate(., change=.[,3]-.[,4]) %>% .[1, "change"]
      
      p1[[7]]$data %>% filter(Crop =="Soy") %>% dplyr::select(starts_with("0.5")) %>% max()
      p1[[7]]$data %>% filter(Crop =="Soy") %>% dplyr::select(starts_with("0.025")| starts_with("0.975")) %>% 
        mutate_at(., .vars=c(1), ~max(.)) %>% mutate_at(., .vars=c(2), ~min(.)) %>% mutate(., change=.[,1]-.[,2]) %>% .[1, "change"]
      p1[[7]]$data %>% filter(Crop =="Soy") %>% dplyr::select(starts_with("0.5")| ID) %>% 
        mutate(top = .[near(ID,1, tol=0.1),1]) %>% mutate(bottom = .[near(ID,-1, tol=0.05),1]) %>% mutate(., change=.[,3]-.[,4]) %>% .[1, "change"]
      
      #SIDI
      p1[[8]]$data %>% filter(Crop =="Winter wheat") %>% dplyr::select(starts_with("0.5")) %>% max()
      p1[[8]]$data %>% filter(Crop =="Winter wheat") %>% dplyr::select(starts_with("0.025")| starts_with("0.975")) %>% 
        mutate_at(., .vars=c(1), ~max(.)) %>% mutate_at(., .vars=c(2), ~min(.)) %>% mutate(., change=.[,1]-.[,2]) %>% .[1, "change"]
      p1[[8]]$data %>% filter(Crop =="Winter wheat") %>% dplyr::select(starts_with("0.5")| ID) %>% 
        mutate(top = .[near(ID,1, tol=0.061),1]) %>% mutate(bottom = .[near(ID,-1, tol=0.1),1]) %>% mutate(., change=.[,3]-.[,4]) %>% .[1, "change"]
      
      p1[[8]]$data %>% filter(Crop =="Corn") %>% dplyr::select(starts_with("0.5")) %>% max()
      p1[[8]]$data %>% filter(Crop =="Corn") %>% dplyr::select(starts_with("0.025")| starts_with("0.975")) %>% 
        mutate_at(., .vars=c(1), ~max(.)) %>% mutate_at(., .vars=c(2), ~min(.)) %>% mutate(., change=.[,1]-.[,2]) %>% .[1, "change"]
      p1[[8]]$data %>% filter(Crop =="Corn") %>% dplyr::select(starts_with("0.5")| ID) %>% 
        mutate(top = .[near(ID,1, tol=0.05),1]) %>% mutate(bottom = .[near(ID,-1, tol=0.2),1]) %>% mutate(., change=.[,3]-.[,4]) %>% .[1, "change"]
      
      p1[[8]]$data %>% filter(Crop =="Soy") %>% dplyr::select(starts_with("0.5")) %>% max()
      p1[[8]]$data %>% filter(Crop =="Soy") %>% dplyr::select(starts_with("0.025")| starts_with("0.975")) %>% 
        mutate_at(., .vars=c(1), ~max(.)) %>% mutate_at(., .vars=c(2), ~min(.)) %>% mutate(., change=.[,1]-.[,2]) %>% .[1, "change"]
      p1[[8]]$data %>% filter(Crop =="Soy") %>% dplyr::select(starts_with("0.5")| ID) %>% 
        mutate(top = .[near(ID,1, tol=0.05),1]) %>% mutate(bottom = .[near(ID,-1, tol=0.1),1]) %>% mutate(., change=.[,3]-.[,4]) %>% .[1, "change"]
      
      #SIEI
      p1[[9]]$data %>% filter(Crop =="Winter wheat") %>% dplyr::select(starts_with("0.5")) %>% max()
      p1[[9]]$data %>% filter(Crop =="Winter wheat") %>% dplyr::select(starts_with("0.025")| starts_with("0.975")) %>% 
        mutate_at(., .vars=c(1), ~max(.)) %>% mutate_at(., .vars=c(2), ~min(.)) %>% mutate(., change=.[,1]-.[,2]) %>% .[1, "change"]
      p1[[9]]$data %>% filter(Crop =="Winter wheat") %>% dplyr::select(starts_with("0.5")| ID) %>% 
        mutate(top = .[near(ID,1, tol=0.061),1]) %>% mutate(bottom = .[near(ID,-1, tol=0.1),1]) %>% mutate(., change=.[,3]-.[,4]) %>% .[1, "change"]
      
      p1[[9]]$data %>% filter(Crop =="Corn") %>% dplyr::select(starts_with("0.5")) %>% max()
      p1[[9]]$data %>% filter(Crop =="Corn") %>% dplyr::select(starts_with("0.025")| starts_with("0.975")) %>% 
        mutate_at(., .vars=c(1), ~max(.)) %>% mutate_at(., .vars=c(2), ~min(.)) %>% mutate(., change=.[,1]-.[,2]) %>% .[1, "change"]
      p1[[9]]$data %>% filter(Crop =="Corn") %>% dplyr::select(starts_with("0.5")| ID) %>% 
        mutate(top = .[near(ID,1, tol=0.05),1]) %>% mutate(bottom = .[near(ID,-1, tol=0.2),1]) %>% mutate(., change=.[,3]-.[,4]) %>% .[1, "change"]
      
      p1[[9]]$data %>% filter(Crop =="Soy") %>% dplyr::select(starts_with("0.5")) %>% max()
      p1[[9]]$data %>% filter(Crop =="Soy") %>% dplyr::select(starts_with("0.025")| starts_with("0.975")) %>% 
        mutate_at(., .vars=c(1), ~max(.)) %>% mutate_at(., .vars=c(2), ~min(.)) %>% mutate(., change=.[,1]-.[,2]) %>% .[1, "change"]
      p1[[9]]$data %>% filter(Crop =="Soy") %>% dplyr::select(starts_with("0.5")| ID) %>% 
        mutate(top = .[near(ID,1, tol=0.07),1]) %>% mutate(bottom = .[near(ID,-1, tol=0.1),1]) %>% mutate(., change=.[,3]-.[,4]) %>% .[1, "change"]
      
      #PNC
      p1[[10]]$data %>% filter(Crop =="Winter wheat") %>% dplyr::select(starts_with("0.5")) %>% max()
      p1[[10]]$data %>% filter(Crop =="Winter wheat") %>% dplyr::select(starts_with("0.025")| starts_with("0.975")) %>% 
        mutate_at(., .vars=c(1), ~max(.)) %>% mutate_at(., .vars=c(2), ~min(.)) %>% mutate(., change=.[,1]-.[,2]) %>% .[1, "change"]
      p1[[10]]$data %>% filter(Crop =="Winter wheat") %>% dplyr::select(starts_with("0.5")| ID) %>% 
        mutate(top = .[near(ID,1, tol=0.061),1]) %>% mutate(bottom = .[near(ID,-1, tol=0.1),1]) %>% mutate(., change=.[,3]-.[,4]) %>% .[1, "change"]
      
      p1[[10]]$data %>% filter(Crop =="Corn") %>% dplyr::select(starts_with("0.5")) %>% max()
      p1[[10]]$data %>% filter(Crop =="Corn") %>% dplyr::select(starts_with("0.025")| starts_with("0.975")) %>% 
        mutate_at(., .vars=c(1), ~max(.)) %>% mutate_at(., .vars=c(2), ~min(.)) %>% mutate(., change=.[,1]-.[,2]) %>% .[1, "change"]
      p1[[10]]$data %>% filter(Crop =="Corn") %>% dplyr::select(starts_with("0.5")| ID) %>% 
        mutate(top = .[near(ID,1, tol=0.05),1]) %>% mutate(bottom = .[near(ID,-1, tol=0.2),1]) %>% mutate(., change=.[,3]-.[,4]) %>% .[1, "change"]
      
      p1[[10]]$data %>% filter(Crop =="Soy") %>% dplyr::select(starts_with("0.5")) %>% max()
      p1[[10]]$data %>% filter(Crop =="Soy") %>% dplyr::select(starts_with("0.025")| starts_with("0.975")) %>% 
        mutate_at(., .vars=c(1), ~max(.)) %>% mutate_at(., .vars=c(2), ~min(.)) %>% mutate(., change=.[,1]-.[,2]) %>% .[1, "change"]
      p1[[10]]$data %>% filter(Crop =="Soy") %>% dplyr::select(starts_with("0.5")| ID) %>% 
        mutate(top = .[near(ID,1, tol=0.07),1]) %>% mutate(bottom = .[near(ID,-1, tol=0.1),1]) %>% mutate(., change=.[,3]-.[,4]) %>% .[1, "change"]
      
      #AREA_MN
      p1[[1]]$data %>% filter(Crop =="Winter wheat") %>% dplyr::select(starts_with("0.025")| starts_with("0.975")| ID)
      p1[[1]]$data %>% filter(Crop =="Winter wheat") %>% dplyr::select(starts_with("0.5")) %>% max()
      p1[[1]]$data %>% filter(Crop =="Winter wheat") %>% dplyr::select(starts_with("0.025")| starts_with("0.975")) %>% 
        mutate_at(., .vars=c(1), ~max(.)) %>% mutate_at(., .vars=c(2), ~min(.)) %>% mutate(., change=.[,1]-.[,2]) %>% .[1, "change"]
      p1[[1]]$data %>% filter(Crop =="Winter wheat") %>% dplyr::select(starts_with("0.5")| ID) %>% 
        mutate(top = .[near(ID,1, tol=0.061),1]) %>% mutate(bottom = .[near(ID,-1, tol=0.05),1]) %>% mutate(., change=.[,3]-.[,4]) %>% .[1, "change"]
      
      p1[[1]]$data %>% filter(Crop =="Corn") %>% dplyr::select(starts_with("0.025")| starts_with("0.975")| ID)
      p1[[1]]$data %>% filter(Crop =="Corn") %>% dplyr::select(starts_with("0.5")) %>% max()
      p1[[1]]$data %>% filter(Crop =="Corn") %>% dplyr::select(starts_with("0.025")| starts_with("0.975")) %>% 
        mutate_at(., .vars=c(1), ~max(.)) %>% mutate_at(., .vars=c(2), ~min(.)) %>% mutate(., change=.[,1]-.[,2]) %>% .[1, "change"]
      p1[[1]]$data %>% filter(Crop =="Corn") %>% dplyr::select(starts_with("0.5")| ID) %>% 
        mutate(top = .[near(ID,1, tol=0.13),1]) %>% mutate(bottom = .[near(ID,-1, tol=0.05),1]) %>% mutate(., change=.[,3]-.[,4]) %>% .[1, "change"]
      
      p1[[1]]$data %>% filter(Crop =="Soy") %>% dplyr::select(starts_with("0.025")| starts_with("0.975")|ID)
      p1[[1]]$data %>% filter(Crop =="Soy") %>% dplyr::select(starts_with("0.5")) %>% max()
      p1[[1]]$data %>% filter(Crop =="Soy") %>% dplyr::select(starts_with("0.025")| starts_with("0.975")) %>% 
        mutate_at(., .vars=c(1), ~max(.)) %>% mutate_at(., .vars=c(2), ~min(.)) %>% mutate(., change=.[,1]-.[,2]) %>% .[1, "change"]
      p1[[1]]$data %>% filter(Crop =="Soy") %>% dplyr::select(starts_with("0.5")| ID) %>% 
        mutate(top = .[near(ID,1, tol=0.11),1]) %>% mutate(bottom = .[near(ID,-1, tol=0.05),1]) %>% mutate(., change=.[,3]-.[,4]) %>% .[1, "change"]
      
      #CONTAG
      p1[[2]]$data %>% filter(Crop =="Winter wheat") %>% dplyr::select(starts_with("0.025")| starts_with("0.975")|ID)
      p1[[2]]$data %>% filter(Crop =="Winter wheat") %>% dplyr::select(starts_with("0.5")) %>% max()
      p1[[2]]$data %>% filter(Crop =="Winter wheat") %>% dplyr::select(starts_with("0.025")| starts_with("0.975")) %>% 
        mutate_at(., .vars=c(1), ~max(.)) %>% mutate_at(., .vars=c(2), ~min(.)) %>% mutate(., change=.[,1]-.[,2]) %>% .[1, "change"]
      p1[[2]]$data %>% filter(Crop =="Winter wheat") %>% dplyr::select(starts_with("0.5")| ID) %>% 
        mutate(top = .[near(ID,1, tol=0.061),1]) %>% mutate(bottom = .[near(ID,-1, tol=0.05),1]) %>% mutate(., change=.[,3]-.[,4]) %>% .[1, "change"]
      
      p1[[2]]$data %>% filter(Crop =="Corn") %>% dplyr::select(starts_with("0.5")) %>% max()
      p1[[2]]$data %>% filter(Crop =="Corn") %>% dplyr::select(starts_with("0.025")| starts_with("0.975")) %>% 
        mutate_at(., .vars=c(1), ~max(.)) %>% mutate_at(., .vars=c(2), ~min(.)) %>% mutate(., change=.[,1]-.[,2]) %>% .[1, "change"]
      p1[[2]]$data %>% filter(Crop =="Corn") %>% dplyr::select(starts_with("0.5")| ID) %>% 
        mutate(top = .[near(ID,1, tol=0.13),1]) %>% mutate(bottom = .[near(ID,-1, tol=0.05),1]) %>% mutate(., change=.[,3]-.[,4]) %>% .[1, "change"]
      
      p1[[2]]$data %>% filter(Crop =="Soy") %>% dplyr::select(starts_with("0.5")) %>% max()
      p1[[2]]$data %>% filter(Crop =="Soy") %>% dplyr::select(starts_with("0.025")| starts_with("0.975")) %>% 
        mutate_at(., .vars=c(1), ~max(.)) %>% mutate_at(., .vars=c(2), ~min(.)) %>% mutate(., change=.[,1]-.[,2]) %>% .[1, "change"]
      p1[[2]]$data %>% filter(Crop =="Soy") %>% dplyr::select(starts_with("0.5")| ID) %>% 
        mutate(top = .[near(ID,1, tol=0.11),1]) %>% mutate(bottom = .[near(ID,-1, tol=0.05),1]) %>% mutate(., change=.[,3]-.[,4]) %>% .[1, "change"]
      
      #ED
      p1[[3]]$data %>% filter(Crop =="Winter wheat") %>% dplyr::select(starts_with("0.5")) %>% max()
      p1[[3]]$data %>% filter(Crop =="Winter wheat") %>% dplyr::select(starts_with("0.025")| starts_with("0.975")) %>% 
        mutate_at(., .vars=c(1), ~max(.)) %>% mutate_at(., .vars=c(2), ~min(.)) %>% mutate(., change=.[,1]-.[,2]) %>% .[1, "change"]
      p1[[3]]$data %>% filter(Crop =="Winter wheat") %>% dplyr::select(starts_with("0.5")| ID) %>% 
        mutate(top = .[near(ID,1, tol=0.1),1]) %>% mutate(bottom = .[near(ID,-1, tol=0.05),1]) %>% mutate(., change=.[,3]-.[,4]) %>% .[1, "change"]
      
      p1[[3]]$data %>% filter(Crop =="Corn") %>% dplyr::select(starts_with("0.5")) %>% max()
      p1[[3]]$data %>% filter(Crop =="Corn") %>% dplyr::select(starts_with("0.025")| starts_with("0.975")) %>% 
        mutate_at(., .vars=c(1), ~max(.)) %>% mutate_at(., .vars=c(2), ~min(.)) %>% mutate(., change=.[,1]-.[,2]) %>% .[1, "change"]
      p1[[3]]$data %>% filter(Crop =="Corn") %>% dplyr::select(starts_with("0.5")| ID) %>% 
        mutate(top = .[near(ID,1, tol=0.13),1]) %>% mutate(bottom = .[near(ID,-1, tol=0.05),1]) %>% mutate(., change=.[,3]-.[,4]) %>% .[1, "change"]
      
      p1[[3]]$data %>% filter(Crop =="Soy") %>% dplyr::select(starts_with("0.5")) %>% max()
      p1[[3]]$data %>% filter(Crop =="Soy") %>% dplyr::select(starts_with("0.025")| starts_with("0.975")) %>% 
        mutate_at(., .vars=c(1), ~max(.)) %>% mutate_at(., .vars=c(2), ~min(.)) %>% mutate(., change=.[,1]-.[,2]) %>% .[1, "change"]
      p1[[3]]$data %>% filter(Crop =="Soy") %>% dplyr::select(starts_with("0.5")| ID) %>% 
        mutate(top = .[near(ID,1, tol=0.11),1]) %>% mutate(bottom = .[near(ID,-1, tol=0.05),1]) %>% mutate(., change=.[,3]-.[,4]) %>% .[1, "change"]
      
      #LPI
      p1[[4]]$data %>% filter(Crop =="Winter wheat") %>% dplyr::select(starts_with("0.025")| starts_with("0.975")|ID)
      p1[[4]]$data %>% filter(Crop =="Winter wheat") %>% dplyr::select(starts_with("0.5")) %>% max()
      p1[[4]]$data %>% filter(Crop =="Winter wheat") %>% dplyr::select(starts_with("0.025")| starts_with("0.975")) %>% 
        mutate_at(., .vars=c(1), ~max(.)) %>% mutate_at(., .vars=c(2), ~min(.)) %>% mutate(., change=.[,1]-.[,2]) %>% .[1, "change"]
      p1[[4]]$data %>% filter(Crop =="Winter wheat") %>% dplyr::select(starts_with("0.5")| ID) %>% 
        mutate(top = .[near(ID,1, tol=0.282),1]) %>% mutate(bottom = .[near(ID,-1, tol=0.1),1]) %>% mutate(., change=.[,3]-.[,4]) %>% .[1, "change"]
      
      p1[[4]]$data %>% filter(Crop =="Corn") %>% dplyr::select(starts_with("0.5")) %>% max()
      p1[[4]]$data %>% filter(Crop =="Corn") %>% dplyr::select(starts_with("0.025")| starts_with("0.975")) %>% 
        mutate_at(., .vars=c(1), ~max(.)) %>% mutate_at(., .vars=c(2), ~min(.)) %>% mutate(., change=.[,1]-.[,2]) %>% .[1, "change"]
      p1[[4]]$data %>% filter(Crop =="Corn") %>% dplyr::select(starts_with("0.5")| ID) %>% 
        mutate(top = .[near(ID,1, tol=0.13),1]) %>% mutate(bottom = .[near(ID,-1, tol=0.1),1]) %>% mutate(., change=.[,3]-.[,4]) %>% .[1, "change"]
      
      p1[[4]]$data %>% filter(Crop =="Soy") %>% dplyr::select(starts_with("0.025")| starts_with("0.975")|ID)
      p1[[4]]$data %>% filter(Crop =="Soy") %>% dplyr::select(starts_with("0.5")) %>% max()
      p1[[4]]$data %>% filter(Crop =="Soy") %>% dplyr::select(starts_with("0.025")| starts_with("0.975")) %>% 
        mutate_at(., .vars=c(1), ~max(.)) %>% mutate_at(., .vars=c(2), ~min(.)) %>% mutate(., change=.[,1]-.[,2]) %>% .[1, "change"]
      p1[[4]]$data %>% filter(Crop =="Soy") %>% dplyr::select(starts_with("0.5")| ID) %>% 
        mutate(top = .[near(ID,1, tol=0.11),1]) %>% mutate(bottom = .[near(ID,-1, tol=0.1),1]) %>% mutate(., change=.[,3]-.[,4]) %>% .[1, "change"]
      
      
      
      