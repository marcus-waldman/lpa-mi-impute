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
    
    # Append in the classification data
    system('rm -r H:\\classification-accuracy')
    system('mkdir H:\\classification-accuracy')
    system('xcopy S:\\classify-accuracy-files\\*.RData H:\\classification-accuracy\\ /J /Y')
    files_vec = list.files("H:/classification-accuracy/", full.names = T)
    classify_df<-
      pblapply(X = 1:500, FUN = function(x){
                                            load(files_vec[x]);
                                            out_x = out_x %>% 
                                              transform(Imputation_Procedure = mapvalues(data_type, from = c("Complete data", "Observed data"), to = c("No Missingness", "FIML")),
                                                        rep = x)
                                            out_x$Imputation_Procedure = as.character(out_x$Imputation_Procedure)
                                            out_x$Imputation_Procedure[out_x$pva==1] = "Amelia"
                                            out_x$Imputation_Procedure[out_x$pva==2] = "PMM"
                                            out_x$Imputation_Procedure[out_x$pva==3] = "CART"
                                            out_x$Imputation_Procedure[out_x$pva==4] = "RF"
                                            return(out_x)
                                } # end function
               )%>%
      rbindlist()

    
    system('mkdir H:\\classification-accuracy\\rand-summary')
    
    
    
    
  
    Processors = 10
    cl<-makeSOCKcluster(Processors)
    doSNOW::registerDoSNOW(cl)
    
#r = 1
    pb <- pbapply::timerProgressBar(max = 500, style = 1, width = getOption("width")/4)
    progress <- function(x){setTimerProgressBar(pb, x)}
    opts <- list(progress = progress)
    tmp_list<-
      foreach(r = 1:500,
              .packages = c("tidyverse","flexclust","data.table"),
              .inorder = TRUE,
              .options.snow = opts) %dopar% {
    
          classify_rep = classify_df %>% filter(rep == r)
          
#x=1
          return_df <- lapply(X = 1:nrow(classify_rep), 
                             FUN = function(x){
                                         hi = classify_rep[x, ] %>% select(starts_with("kappa")) %>% as.integer() %>% matrix(nrow = 3)
                                         if(sum(is.na(hi))==0){
                          
                                           Rand = randIndex(hi %>% as.table(), correct = FALSE)
                          
                                         } else {
                                           Rand = NA
                                         }
                                         return(with(classify_rep,
                                                    data.frame(Imputation_Procedure = Imputation_Procedure[x],
                                                               rep = rep[x],
                                                               z = z[x],
                                                               pm = pm[x],
                                                               m = m[x],
                                                               pva = pva[x],
                                                               Rand = Rand)))
                                   }
                             ) %>% rbindlist()
           
                  save(return_df, file = paste0("H:/classification-accuracy/rand-summary/rand-summary-rep",r,".RData"))

                  return(NULL)
          
              } #end foreach
    stopCluster(cl)
    
    
    
    # Copy over
    system('xcopy H:\\classification-accuracy\\rand-summary\\*.RData S:\\classify-accuracy-files\\rand-summary /J /Y')
    
