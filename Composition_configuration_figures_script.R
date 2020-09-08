library(pacman)

p_load(tidyverse, INLA, sf, spdep, rgeos, knitr, ggpubr,doParallel, foreach,tidyselect, viridis)


wd<-getwd()

source(paste0(wd,'/Scripts/SA_functions.R'))

##################################################
### Create summary plots for landscape metrics ###
##################################################

  ivars<-c("LSM_AREA_MN_ALL", "LSM_CONTAG_ALL","LSM_ED_ALL","LSM_LPI_ALL",
           "LSM_RICH_ALL","LSM_SHDI_ALL","LSM_SHEI_ALL","LSM_SIDI_ALL","LSM_SIEI_ALL", "PNC")
  
  metrics<-ivars
  crops<-c("corn","soy","wwheat")
  y.min<-(-0.15)
  y.max<-(0.2)
  x.max<-c(3.25,3.25,3.25,3.25,3.25,3.25,3.25,3.25,3.25,3.25)
  x.min<-c(-3.25,-3.25,-3.25,-3.25,-3.25,-3.25,-3.25,-3.25,-3.25,-3.25)
  
  p1<-plot.lm.summ(metrics, crops, y.min, y.max, x.min, x.max)
  
  
    #build a 5 x 2 plot by combining the output for all landscape metrics and crops
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
    ggsave("composition.png", figure1, width = 5, height = 6)
    dev.off()
    
    #build a 2 x 2 plot by combining the output for all configuration metrics and crops
    figure2<-ggarrange(p1[[1]], p1[[4]] + theme( axis.title.y = element_blank()), p1[[2]] , p1[[3]] + theme( axis.title.y = element_blank()),
                       ncol=2,nrow=2, common.legend = TRUE, legend = "bottom", 
                       labels="auto",hjust = c(-3, -1,-3,-1), vjust =c(11,11,11,11) )
    annotate_figure(figure2)
    ggsave("configuration.png", figure2, width = 5, height = 4)
    dev.off()

########################################################################################
### Create simplified climate curve summary plot across all metric-crop combinations ###
########################################################################################
    
  metrics<-ivars
  crops<-c("corn","soy","wwheat")
  y.min<- c(-0.2, -0.2, -0.55, -0.1)
  y.max<-c(0.1, 0.3, 0.3, 0.3)
  x.max<- c(3.5, 3.5, 3.5, 3.5)
  x.min<- c(-2.5, -2.5, -2.5, -2.5)
  
  p3<-plot.cl.summ.simple(metrics, crops, y.min, y.max, x.min, x.max)
  
    #build a 2x2 plot of climate curves
    figure3<-ggarrange(p3[[2]], p3[[3]] + theme( axis.title.y = element_blank()), 
                       p3[[1]] , p3[[4]] + theme( axis.title.y = element_blank()),
                       ncol=2,nrow=2, common.legend = TRUE, legend = "bottom", 
                       labels="auto",hjust = c(-3, -1,-3,-1), vjust =c(11,11,11,11) )
    annotate_figure(figure3)
    ggsave("climate_summary.png", figure3, width = 5, height = 4)


#########################################################################
### Create plot of climate curves across all metric-crop combinations ###
#########################################################################
    
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
    

###############################################
### Create a summary interaction model plot ###
###############################################
    
    metrics<-c('SHDIxED','SHDIxMNAREA','RICHxED','RICHxMNAREA')
    crops<-c('corn','soy','wwheat')
    y.min<-(-0.2)
    y.max<-(0.325)
    x.max<-c(3.25)
    x.min<-c(-3.25)
    
    p<-plot.interaction(metrics, crops, y.min, y.max, x.min, x.max)
    
    #build a 3x1 plot of interaction curves
    
    
    figure4<-ggarrange(p[[1]] + theme( plot.title=element_blank()), 
                       p[[3]] + theme( axis.title.y = element_blank(), plot.title=element_blank()), 
                       p[[2]] + theme( axis.title.y = element_blank(),  plot.title=element_blank()), 
                       ncol=3,nrow=1, common.legend = TRUE, legend = "top", 
                       labels="auto",hjust = c(-2.5, -1,-1), vjust =c(17,17,17) )
    annotate_figure(figure4)
    ggsave("interactions.png", figure4, width = 6.5, height = 3)
    