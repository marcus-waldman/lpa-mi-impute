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
    wd = "C:/Users/marcu/Documents/lpa-mi-stage0"
    #datfiles_wd = "D:/datfiles"
    datfiles_wd = paste0(wd,"/mi-impute-stage1-create Jan 11 2019 19 41")
    outfiles_wd = paste0(wd,"/stage1")
    
    # Number of classes fit 
    K_vec = 3;
    tmax = 1/6; #Maximum (expected) number of minutes to fit a single model before stopping.
    
    # Set the working directory
    setwd(wd)
    
    # Load the working environment
    load(paste0(wd,"/Rdata-mi-impute-stage1-create Jan 11 2019 19 41/environment.RData"))
    
    

    outwd = outfiles_wd # Where wouldyou like to save the out files
    datwd = datfiles_wd #Where is the home directory of the dat files
    data_type = "Imputed data"
    rep_vec = 1:500
    z_vec = 1
    data_conditions
    pm_vec = 1
    pva_vec = 3
    M = 20
    dk_vec = 0
    starts0 = 60
    multicore = FALSE

    # tracker_df = gen_tracker_and_outdirs(outwd = outfiles_wd,
    #                                       datwd = datfiles_wd,
    #                                       rdatawd = "H:/Google Drive/diss-drive",
    #                                       data_type = "Imputed data",
    #                                       rep_vec = 1:500,
    #                                       z_vec = 1,
    #                                       pm_vec = 1,
    #                                       pva_vec = methods_list$pva_vec,
    #                                       data_conditions = transform(data_conditions, z = 1:nrow(data_conditions)),
    #                                       M = 20)
    # 
    # save(tracker_df, file = "tracker.RData")
    # 
    # 
    # tracker_df = tracker2_fit_til_convergence(data_conditions = data_conditions,
    #                                         tracker_df = tracker_df,
    #                                         Processors = Processors,
    #                                         tmax = tmax)
    # 
    # 
    # save(tracker_df, file = "tracker.RData")
# 
# 
# 
#     tracker_df$converged = FALSE
#     tracker_df$converged[tracker_df$normal==TRUE & tracker_df$poor_Rcond==FALSE] = TRUE
#     tracker_df = tracker2_switch(tracker_df = tracker_df,
#                              data_conditions = data_conditions)
# 
#     save(tracker_df, file = "tracker.RData")

    load("tracker.RData")
    
    
    
    tracker_df = tracker2_save_readMplus(tracker_df=tracker_df, estwd = "C:/Users/marcu/Documents/lpa-mi-stage0")
    save(tracker_df, file = "tracker.RData")
    
    
   

      
  