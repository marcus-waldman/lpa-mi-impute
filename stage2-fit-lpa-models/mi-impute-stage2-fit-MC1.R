#  .rs.restartR()
    
    rm(list = ls())
    
    library(plyr)
    library(tidyverse)
    library(lpa.mi.src)
    library(snow)
    library(doSNOW)
    library(foreach)
    library(data.table)
    
    computer_name = "MC1"
    wd = "H:"
    dat_wd = paste0(wd,"/dat-files")
    out_wd = paste0(wd,"/out-files")
    rdata_wd = paste0(wd,"/Google Drive/diss-drive")
    pingpong_wd = paste0("D:/Dropbox/Dissertation/ping-pong")
    
    Processors = 4
    cl<-makeSOCKcluster(Processors)
    doSNOW::registerDoSNOW(cl)    
    
    # Number of classes fit 
    K_vec = 3;

    # Set the working directory
    setwd(wd)
    
    # Load the working environment
    load("environment-mi-impute-stage0 Jan 06 2019 13 29.RData")
    #rm(methods_list)
    

    z_vec = 1:40
    pm_vec = 1
    pva_vec = 1:4
    M = 100
    dk_vec = 0
    starts0 = 20

repseq = seq(500,1,by=-1)

for (rep in repseq){
          
          # Check if should break, or serve the ping.
      if( !(paste0("rep",rep,".csv")%in%list.files(path = pingpong_wd)) ){
          
          write.csv(x = data.frame(computer = computer_name, finished = FALSE), file = paste0(pingpong_wd,"/rep",rep,".csv"), row.names = FALSE)    
          complete_tracker = gen_tracker_and_outdirs(outwd = out_wd,
                                                    datwd = dat_wd,
                                                    rdatawd = rdata_wd,
                                                    data_type = "Complete data",
                                                    rep_vec = rep,
                                                    z_vec = z_vec,
                                                    pm_vec = pm_vec,
                                                    pva_vec = pva_vec,
                                                    data_conditions = transform(data_conditions, z = 1:nrow(data_conditions)),
                                                    M = M, 
                                                    starts0 = starts0)
         
          
          
          obs_tracker = gen_tracker_and_outdirs(outwd = out_wd,
                                                datwd = dat_wd,
                                                rdatawd = rdata_wd,
                                                data_type = "Observed data",
                                                rep_vec = rep,
                                                z_vec = z_vec,
                                                pm_vec = pm_vec,
                                                pva_vec = pva_vec,
                                                data_conditions = transform(data_conditions, z = 1:nrow(data_conditions)),
                                                M = M, 
                                                starts0 = starts0)
          
          imp_tracker = gen_tracker_and_outdirs(outwd = out_wd,
                                                     datwd = dat_wd,
                                                     rdatawd = rdata_wd,
                                                     data_type = "Imputed data",
                                                     rep_vec = rep,
                                                     z_vec = z_vec,
                                                     pm_vec = pm_vec,
                                                     pva_vec = pva_vec,
                                                     data_conditions = transform(data_conditions, z = 1:nrow(data_conditions)),
                                                     M = M, 
                                                     starts0 = starts0)
          
          tracker_df = rbind(complete_tracker, obs_tracker, imp_tracker)
          tracker_df$nn = 1:nrow(tracker_df)
      
          if(length(unique(tracker_df$nn))!=length(tracker_df$nn)){stop("nn identifier in tracker_df is not unique")}
          #save(tracker_df, file = paste0("tracker-rep", rep, ".RData"))
      
      
          tracker_df = tracker2_fit_til_convergence(data_conditions = data_conditions,
                                                  tracker_df = tracker_df,
                                                  cl = cl)
  
          
          if(length(unique(tracker_df$nn))!=length(tracker_df$nn)){stop("nn identifier in tracker_df is not unique")}
          #save(tracker_df, file = paste0("tracker-rep", rep, ".RData"))
      
          
          tracker_df$converged = FALSE
          tracker_df$converged[tracker_df$normal==TRUE & tracker_df$poor_Rcond==FALSE] = TRUE
          tracker_df = tracker2_switch(tracker_df = tracker_df,
                                   data_conditions = data_conditions, 
                                   cl = cl)

                
          if(length(unique(tracker_df$nn))!=length(tracker_df$nn)){stop("nn identifier in tracker_df is not unique")}
          #save(tracker_df, file = paste0("tracker-rep", rep, ".RData"))
          
          tracker_df = tracker2_save_readMplus(tracker_df=tracker_df, estwd = wd, cl = cl)
          
          if(length(unique(tracker_df$nn))!=length(tracker_df$nn)){stop("nn identifier in tracker_df is not unique")}
          
          
          save(tracker_df, file = paste0("tracker-rep", rep, ".RData"))
          write.csv(x = data.frame(computer = computer_name, finished = TRUE), file = paste0(pingpong_wd,"/rep",rep,".csv"), row.names = FALSE)    
          
          
  }# end if
} #rep in repseq

stopCluster(cl)

    
   

      
  