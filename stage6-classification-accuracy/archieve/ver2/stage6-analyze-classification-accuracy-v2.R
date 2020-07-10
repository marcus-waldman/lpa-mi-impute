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
    system('mkdir H:\\classify-accuracy-files\\rand-summary')
    system('xcopy S:\\classify-accuracy-files\\rand-summary\\*.RData H:\\classify-accuracy-files\\rand-summary /J /Y')
    
    files_vec = list.files("H:/classify-accuracy-files/rand-summary", full.names = T)
    Rand_df<-
      pblapply(X = 1:500, 
               FUN = function(x){
                  load(files_vec[x]);
                  summary_df = return_df %>% group_by(Imputation_Procedure,rep,z) %>% summarise(Rand = mean(Rand))
                  return(summary_df)
               }) %>% 
      rbindlist()

            
    
    
    setwd(environment_wd)
    load(file = "labels_for_conditions.Rdata")
    load(file = "data-conditions.Rdata")
    
    
    # Load in the data
    setwd(results_wd)
    load(file ="summaries-combined-results-lpa-mi-impute.RData")
    
    
    summaries_combined_df<- summaries_combined_df %>% 
                            left_join(Rand_df, by = c("rep","z","Imputation_Procedure"))
    
    agg_df = summaries_combined_df %>% 
      left_join(conditions_labels, by = "z") %>% 
      group_by(Imputation_Procedure,Mixing,Sample.Size,Separation) %>% 
      summarise(Rand = mean(Rand, na.rm = T)) %>% 
      filter(Imputation_Procedure != "PMM" & Imputation_Procedure != "No Missingness")
    
    ggplot(agg_df %>% filter(startsWith(as.character(Sample.Size), "Large")), aes(x = Imputation_Procedure, y = Rand)) + geom_point() + facet_grid(Separation~Mixing)
    
    ggplot(agg_df %>% filter(startsWith(as.character(Sample.Size), "Small")), aes(x = Imputation_Procedure, y = Rand)) + geom_point() + facet_grid(Separation~Mixing)
    
    # Same effect as before. I think we have to note that the mixtures are getting more separated due to nonresponse bias and so that is why FIML is outperforming the others.