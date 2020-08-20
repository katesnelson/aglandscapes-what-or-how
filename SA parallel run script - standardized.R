library(ggplot2)
library(dplyr)
library(sf)
library(INLA)
library(raster)
library(spdep)
library(rgeos)
library(knitr)
library(ggpubr)
library(doParallel)
library(foreach)

wd<-getwd()
source(paste0(wd,'/Scripts/SA_functions.R'))

    # d<-readRDS(paste0(wd,'/Data/cornstd.RDS'))
    # d<-readRDS(paste0(wd,'/Data/soystd.RDS'))
    # d<-readRDS(paste0(wd,'/Data/wwheatstd.RDS'))

ivars<-c("LSM_AREA_MN_ALL", "LSM_CONTAG_ALL","LSM_ED_ALL","LSM_LPI_ALL","LSM_RICH_ALL","LSM_SHDI_ALL","LSM_SHEI_ALL","LSM_SIDI_ALL","LSM_SIEI_ALL", "PNC")
datafiles<-c("corn_panel", "soy_panel","wwheat_panel")
crops<-c("corn","soy","wwheat")

foreach (j=1:length(datafiles))%do%{
prep.data.std(file=datafiles[j], savename=paste0(crops[j],"std"), crop=crops[j], projection="+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m no_defs") #102003 epsg code for USA albers equal area conic
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



ivars<-c("LSM_AREA_MN_ALL", "LSM_CONTAG_ALL","LSM_ED_ALL","LSM_LPI_ALL","LSM_RICH_ALL","LSM_SHDI_ALL","LSM_SHEI_ALL","LSM_SIDI_ALL","LSM_SIEI_ALL", "PNC")

metrics<-ivars
crops<-c("corn","soy","wwheat")
y.min<-(-0.15)
y.max<-(0.2)
x.max<-c(3.25,3.25,3.25,3.25,3.25,3.25,3.25,3.25,3.25,3.25)
x.min<-c(-3.25,-3.25,-3.25,-3.25,-3.25,-3.25,-3.25,-3.25,-3.25,-3.25)

p1<-plot.lm.summ(metrics, crops, y.min, y.max, x.min, x.max)


#build a 5 x 2 plot by combining the output for all metrics and crops
figure1<-ggarrange(p1[[1]], p1[[2]] + theme( axis.title.y = element_blank()), p1[[3]] , p1[[4]] + theme( axis.title.y = element_blank()),
                   p1[[5]] , p1[[6]] + theme( axis.title.y = element_blank()), p1[[7]] , p1[[8]] + theme( axis.title.y = element_blank()),
                   p1[[9]] , p1[[10]] + theme( axis.title.y = element_blank()),
                   ncol=2,nrow=5, common.legend = TRUE, legend = "bottom")
annotate_figure(figure1)
dev.off()

#build a 3x2 plot by combining the output for all composition metrics and crops
figure1<-ggarrange(p1[[5]] , p1[[6]] + theme( axis.title.y = element_blank()), 
                   p1[[7]] , p1[[8]] + theme( axis.title.y = element_blank()),
                   p1[[9]] , p1[[10]] + theme( axis.title.y = element_blank()),
                   ncol=2,nrow=3, common.legend = TRUE, legend = "bottom", 
                   labels="auto",  hjust = c(-3, -1,-3,-1,-3,-3), vjust =c(11,11,11,11,11,11)) #c("a","b","c","d","e","f")
annotate_figure(figure1)
ggsave("composition_new.png", figure1, width = 5, height = 6)
dev.off()

#build a 2 x 2 plot by combining the output for all configuration metrics and crops
figure2<-ggarrange(p1[[1]], p1[[2]] + theme( axis.title.y = element_blank()), p1[[3]] , p1[[4]] + theme( axis.title.y = element_blank()),
                   ncol=2,nrow=2, common.legend = TRUE, legend = "bottom", 
                   labels="auto",hjust = c(-3, -1,-3,-1), vjust =c(11,11,11,11) )
annotate_figure(figure2)
ggsave("configuration_new.png", figure2, width = 5, height = 4)
dev.off()


##AND STANDARDIZED PLOTS
    # metrics<-ivars
    # crops<-c("corn","soy","wwheat")
    # y.min<-(-0.15)
    # y.max<-(0.2)
    # x.max<-c(3,3,3,3,3,3,3,3,3,3)
    # x.min<-c(-3,-3,-3,-3,-3,-3,-3,-3,-3,-3)
    # 
    # p1<-plot.lm.summ.std(metrics, crops, y.min, y.max, x.min, x.max)
    # 
    # 
    # figure1<-ggarrange(p1[[1]], p1[[2]] + theme( axis.title.y = element_blank()), p1[[3]] , p1[[4]] + theme( axis.title.y = element_blank()),
    #                    p1[[5]] , p1[[6]] + theme( axis.title.y = element_blank()), p1[[7]] , p1[[8]] + theme( axis.title.y = element_blank()),
    #                    p1[[9]] , p1[[10]] + theme( axis.title.y = element_blank()),
    #                    ncol=2,nrow=5, common.legend = TRUE, legend = "bottom")
    # annotate_figure(figure1)

#plot simplified climate curve summary across all metric-crop combinations
  metrics<-ivars
  crops<-c("corn","soy","wwheat")
  y.min<- c(-0.2, -0.2, -0.55, -0.1)
  y.max<-c(0.1, 0.3, 0.3, 0.3)
  x.max<- c(3.5, 3.5, 3.5, 3.5)
  x.min<- c(-2.5, -2.5, -2.5, -2.5)
  
  p3<-plot.cl.summ.simple(metrics, crops, y.min, y.max, x.min, x.max)
  
  #build a 2x2 plot of climate curves
  figure3<-ggarrange(p3[[1]], p3[[4]] + theme( axis.title.y = element_blank()), 
                     p3[[2]] , p3[[3]] + theme( axis.title.y = element_blank()),
                    ncol=2,nrow=2, common.legend = TRUE, legend = "bottom", 
                    labels="auto",hjust = c(-3, -1,-3,-1), vjust =c(11,11,11,11) )
  annotate_figure(figure3)
  ggsave("climate_new.png", figure3, width = 5, height = 4)
  

#plot climate curves across all metric-crop combinations
metrics<-ivars
crops<-c("corn","soy","wwheat")
y.min<-data.frame("corn" = c(-0.1, -0.25, -0.6, -0.1), "soy" = c( -0.2, -0.15, -0.6, -0.1), "wwheat"= c(-0.4, -0.4, -0.2, -0.1))
y.max<-data.frame("corn"= c(0.1, 0.5, 0.35, 0.5), "soy"= c(0.1, 0.1, 0.2, 0.25), "wwheat"= c(0.1, 00.4, 0.2, 0.25))
x.max<-data.frame("corn" = c(1750, 5000, 650, 1),"soy" = c(1750, 5000, 650, 1),"wwheat" = c(2000, 8000, 650, 1))
x.min<-data.frame("corn"= c(0, 1500, 0, 0),"soy"= c(0, 1500, 0, 0),"wwheat"= c(0, 3500, 0, 0))
                  
p3<-plot.lm.summ.cl(metrics, crops, y.min, y.max, x.min, x.max)

#build a 3x4 plot of climate curves
figure1<-ggarrange(p3[[1]][[1]], p3[[1]][[2]]+ theme( axis.title.y = element_blank()), p3[[1]][[3]]+ theme( axis.title.y = element_blank()), p3[[1]][[4]]+ theme( axis.title.y = element_blank()),
                   p3[[2]][[1]], p3[[2]][[2]]+ theme( axis.title.y = element_blank()), p3[[2]][[3]]+ theme( axis.title.y = element_blank()), p3[[2]][[4]]+ theme( axis.title.y = element_blank()),
                   p3[[3]][[1]], p3[[3]][[2]]+ theme( axis.title.y = element_blank()), p3[[3]][[3]]+ theme( axis.title.y = element_blank()), p3[[3]][[4]]+ theme( axis.title.y = element_blank()),
                   ncol=4,nrow=3, common.legend = TRUE, legend = "bottom")
annotate_figure(figure1)

#AND STANDARDIZED
    # metrics<-ivars
    # crops<-c("corn","soy","wwheat")
    # y.min<-data.frame("corn" = c(-0.1, -0.25, -0.6, -0.1), "soy" = c( -0.2, -0.15, -0.6, -0.1), "wwheat"= c(-0.4, -0.4, -0.2, -0.1))
    # y.max<-data.frame("corn"= c(0.1, 0.5, 0.35, 0.5), "soy"= c(0.1, 0.1, 0.2, 0.25), "wwheat"= c(0.1, 00.4, 0.2, 0.25))
    # x.max<-data.frame("corn" = c(1750, 5000, 650, 1),"soy" = c(1750, 5000, 650, 1),"wwheat" = c(2000, 8000, 650, 1))
    # x.min<-data.frame("corn"= c(0, 1500, 0, 0),"soy"= c(0, 1500, 0, 0),"wwheat"= c(0, 3500, 0, 0))
    # 
    # p3<-plot.lm.summ.cl.std(metrics, crops, y.min, y.max, x.min, x.max)
    # 
    # #build a 3x4 plot of climate curves
    # figure1<-ggarrange(p3[[1]][[1]], p3[[1]][[2]]+ theme( axis.title.y = element_blank()), p3[[1]][[3]]+ theme( axis.title.y = element_blank()), p3[[1]][[4]]+ theme( axis.title.y = element_blank()),
    #                    p3[[2]][[1]], p3[[2]][[2]]+ theme( axis.title.y = element_blank()), p3[[2]][[3]]+ theme( axis.title.y = element_blank()), p3[[2]][[4]]+ theme( axis.title.y = element_blank()),
    #                    p3[[3]][[1]], p3[[3]][[2]]+ theme( axis.title.y = element_blank()), p3[[3]][[3]]+ theme( axis.title.y = element_blank()), p3[[3]][[4]]+ theme( axis.title.y = element_blank()),
    #                    ncol=4,nrow=3, common.legend = TRUE, legend = "bottom")
    # annotate_figure(figure1)


#get model fit statistics
model.results.summ(data="cornfinal", crop="corn", metrics)
model.results.summ(data="soyfinal", crop="soy", metrics)
model.results.summ(data="wwheatfinal", crop="wwheat", metrics)


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
      
      
#Get max impact (median effect), total change (conservative, using confidence bounds), and -1 to +1 sd change (median effect)
      
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
      
      
      #Compute R2 as variance explained (sum of squares of full model in comparison to null model) http://www.stat.columbia.edu/~gelman/research/published/rsquared.pdf, http://www.stat.columbia.edu/~gelman/research/unpublished/bayes_R2.pdf
      
      
      
      