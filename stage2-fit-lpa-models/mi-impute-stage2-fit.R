#  .rs.restartR()
    
    rm(list = ls())
    
    library(plyr)
    library(tidyverse)
    library(lpa.mi.src)
    library(doParallel)
    library(foreach)
    library(doRNG)
    library(pbapply)
    library(data.table)
    
    Processors = 15
    wd = "H:"
    dat_wd = paste0(wd,"/dat-files")
    out_wd = paste0(wd,"/out-files")
    rdata_wd = paste0(wd,"/Google Drive/diss-drive")
    pingpong_wd = paste0("D:/Dropbox/Dissertation/ping-pong")
    
    # Number of classes fit 
    K_vec = 3;
    tmax = 1/6; #Maximum (expected) number of minutes to fit a single model before stopping.
    
    # Set the working directory
    setwd(wd)
    
    # Load the working environment
    load(paste0(wd,"/rdata-files/environment - mi-impute-stage1-create Jan 06 2019 22 50.RData"))
    rm(methods_list)
    

    rep_vec = 1
    z_vec = 1:Ztot
    pm_vec = 1
    pva_vec = c(1:4)
    M = 100
    dk_vec = 0
    starts0 = 30
    multicore = FALSE

    
    for (rep in 1:500){
          
          # Check if should break, or serve the ping.
          if(paste0("rep",rep,".csv")%in%list.files(path = pingpong_wd)){
            paste0("<<>> Ping has a pong at rep ", rep,". Ending. <<>>")
            break
          }
          
          write.csv(x = rep, file = paste0("rep",rep,".csv"))    
      
          complete_tracker = gen_tracker_and_outdirs(outwd = out_wd,
                                                    datwd = dat_wd,
                                                    rdatawd = rdata_wd,
                                                    data_type = "Complete data",
                                                    rep_vec = rep_vec,
                                                    z_vec = z_vec,
                                                    pm_vec = pm_vec,
                                                    pva_vec = pva_vec,
                                                    data_conditions = transform(data_conditions, z = 1:nrow(data_conditions)),
                                                    M = M)
          
          obs_tracker = gen_tracker_and_outdirs(outwd = out_wd,
                                                datwd = dat_wd,
                                                rdatawd = rdata_wd,
                                                data_type = "Observed data",
                                                rep_vec = rep_vec,
                                                z_vec = z_vec,
                                                pm_vec = pm_vec,
                                                pva_vec = pva_vec,
                                                data_conditions = transform(data_conditions, z = 1:nrow(data_conditions)),
                                                M = M)
          
          imp_tracker = gen_tracker_and_outdirs(outwd = out_wd,
                                                     datwd = dat_wd,
                                                     rdatawd = rdata_wd,
                                                     data_type = "Imputed data",
                                                     rep_vec = rep_vec,
                                                     z_vec = z_vec,
                                                     pm_vec = pm_vec,
                                                     pva_vec = pva_vec,
                                                     data_conditions = transform(data_conditions, z = 1:nrow(data_conditions)),
                                                     M = M)
          
          tracker_df = rbind(complete_tracker, obs_tracker, imp_tracker)
      
          save(tracker_df, file = paste0("tracker-rep", rep, ".RData"))
      
      
          tracker_df = tracker2_fit_til_convergence(data_conditions = data_conditions,
                                                  tracker_df = tracker_df,
                                                  Processors = Processors,
                                                  tmax = tmax)
      
      
          save(tracker_df, file = paste0("tracker-rep", rep, ".RData"))
      
      
      
          tracker_df$converged = FALSE
          tracker_df$converged[tracker_df$normal==TRUE & tracker_df$poor_Rcond==FALSE] = TRUE
          tracker_df = tracker2_switch(tracker_df = tracker_df,
                                   data_conditions = data_conditions)
      
          save(tracker_df, file = paste0("tracker-rep", rep, ".RData"))
          
          tracker_df = tracker2_save_readMplus(tracker_df=tracker_df, estwd = "C:/Users/marcu/Documents/lpa-mi-stage0")
          save(tracker_df, file = paste0("tracker-rep", rep, ".RData"))
    } #rep in 1:500
    
    
   

      
  