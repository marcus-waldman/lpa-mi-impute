    rm(list = ls())
    gc()
    
    library(plyr)
    library(tidyverse)
    library(data.table)
    library(pbapply)
    library(stringr)
    library(lpa.mi.src)
    library(tangram)
    library(dplyr)
    library(qwraps2)
    library(Hmisc)
    library(flexclust)
    
    library(doParallel)
    library(foreach)
    library(doRNG)
    require(snow)
    require(doSNOW)
    require(foreach)
    require(pbapply)
    
    
    
    # Directories
    dropbox_wd = "D:/Dropbox"
    #dropbox_wd = "C:/Users/marcu/Dropbox"
    results_wd = paste0(dropbox_wd, "/Dissertation/lpa-mi-impute/stage4c-combine-results")
    stage6_wd = paste0(dropbox_wd, "/Dissertation/lpa-mi-impute/stage6-classification-accuracy")
    environment_wd = paste0(dropbox_wd,"/Dissertation/environmental-variables/")
    
    
    setwd(stage6_wd)
    
    system('rm -r H:\\classify-accuracy-files')
    system('mkdir H:\\classify-accuracy-files')
    system('xcopy S:\\classify-accuracy-files\\*.Rdata H:\\classify-accuracy-files /J /Y')
    
    
    classify_df<-
      pblapply(X = 1:500, 
               FUN = function(x){
                 load(paste0("H:/classify-accuracy-files/v3-classify-accuracy-rep",x,".RData"));
                 out_x = out_x %>% 
                   transform(Imputation_Procedure = mapvalues(data_type, from = c("Complete data", "Observed data"), 
                                                                          to = c("No Missingness", "FIML")),
                             rep = x)
                 out_x$Imputation_Procedure = as.character(out_x$Imputation_Procedure)
                 out_x$Imputation_Procedure[out_x$pva==1] = "Amelia"
                 out_x$Imputation_Procedure[out_x$pva==2] = "PMM"
                 out_x$Imputation_Procedure[out_x$pva==3] = "CART"
                 out_x$Imputation_Procedure[out_x$pva==4] = "RF"
                  return(out_x)
               }) %>% rbindlist()
  
            
    
    
    setwd(environment_wd)
    load(file = "labels_for_conditions.Rdata")
    load(file = "data-conditions.Rdata")
    
    
    # Load in the data
    setwd(results_wd)
    load(file ="summaries-combined-results-lpa-mi-impute.RData")
    
    
    classify_df<- classify_df %>% 
      left_join(summaries_combined_df %>% select(rep,z,Imputation_Procedure,Q,MM,ARIV,ARIU,KL,se_KL,points_montecarlo), 
                by = c("rep","z","Imputation_Procedure")) %>% 
      left_join(conditions_labels, by = "z")
    
    classify_df$keep = 1
    classify_df$keep[classify_df$summary_type=="EAPP"]=0
    classify_df$keep[classify_df$summary_type=="pooled, log odds ratio scaled"]=0
    classify_df$keep[classify_df$pva==2]=0
    
  
    agg_df = classify_df %>% 
      filter(keep == 1) %>% 
      group_by(Imputation_Procedure,summary_type, Mixing,Sample.Size,Separation) %>% 
      summarise(kl_cprob = mean(kl_cprob, na.rm = T)) 
    
    ggplot(agg_df %>% filter(startsWith(as.character(Sample.Size), "Large")  ), aes(x = Imputation_Procedure, y = -1*kl_cprob)) + geom_point() + facet_grid(Separation~Mixing)
    
    ggplot(agg_df %>% filter(startsWith(as.character(Sample.Size), "Small")), aes(x = Imputation_Procedure, y = -1*kl_cprob)) + geom_point() + facet_grid(Separation~Mixing)
    