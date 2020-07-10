    rm(list = ls())
    
    library(plyr)
    library(tidyverse)
    library(lpa.mi.src)
    library(doParallel)
    library(foreach)
    library(doRNG)
    
    Processors = 50
    Replications = 10
    #KKmax = 4;
    
    # Acer-laptop
      #dropbox_wd =  "C:/Users/marcu/Dropbox"  #  Dropbox folder on local machine
      #temp_dir = "C:/Users/marcu/Documents/lpa-mi-impute" #  Temporary directory for local maching (will automaticall add new temporary folder)
    # Queen Mary's Revenge
      #dropbox_wd =  "D:/Dropbox"  #  Dropbox folder on local machine
      #temp_dir = paste0(LETTERS[seq(26-Processors+1,26)],":") #  Temporary directory for local maching (will automaticall add new temporary folder)        
    # Xi-GSU
      #dropbox_wd = "C:/Users/mwaldman1/Dropbox"
      #temp_dir = "C:/Users/mwaldman1/Documents/lpa-mi-impute"
    # Big Bertha
      #dropbox_wd =  "D:/Dropbox/Dropbox"  #  Dropbox folder on local machine
      #temp_dir = "D:/lpa-mi-impute" #  Temporary directory for local maching (will automaticall add new temporary folder)        
    # Gaby
      dropbox_wd = "Z:/Dropbox"
      temp_dir = c(rep("W:/", ceiling(Processors/3)), 
                 rep("X:/", ceiling(Processors/3)), 
                 rep("Y:/", Processors - 2*ceiling(Processors/3)))
    
    # Nuke directory
      lpa.mi.src::nuke_dirs(temp_dir)
        
    # Create folder to store results in main directory
      setwd(paste(dropbox_wd, "/Dissertation/results", sep = ""))
      temp_fname = paste("sim-lpa-mi-impute", format(Sys.time(), "%b %d %Y %H %M"), sep = " ")
      dir.create(temp_fname)
      setwd(temp_fname)
      results_wd = getwd()
        
    
    # Fully-crossed simulation conditions
    pctmiss_vec = c(0.5)     #  Percent of observations with *at least* one missing value (indicator or missing data correlate) -> if other specification needed one must modify list_inds.miss generation in get_obs_data.R
    N_vec = c(2E3)                #  Sample sizes within each replication
    J_Y_vec = c(3)                #  Number of latent class indicators
    J_Xcom_vec = c(1)             #  Number of complete data missing data correlates
    J_Xinc_vec = c(0)             #  Number of missing data correlates which themselves contain missing data
    MD_vec = c(10/3)               #  WARNING! MD = 1 RESULTS IN ERROR! Mahalanobis distances between classes
    pi_list = list(               #  Number of classes, K (maximum of 3)
      A = rep(1/3,3)#, 
      #B = c(0.425,0.425,0.15) #, 
      #C = rep(1/4,4)
    )
    rho_YX_vec = c(0.4)          #  Strength of missing data correlates
    C_modifies_YX_vec = c(TRUE)  # Whether or not the Y-X relationships are different between latent classes 
    t_rotate_vec = c(0,pi/4)     # Angle for rotation matrix for LPA indicators 
    
    # Define the imputation alogrithms to compare
    methods_list = list(procedure = list(NULL), name = list(NULL), args = list(NULL))
    methods_list$procedure[[1]] = "stratamelia";    methods_list$name[[1]] = "mg-Amelia";    methods_list$args[[1]] = list(m = 1, p2s = 0, tolerance = 1E-4)
    methods_list$procedure[[2]] = "mice";           methods_list$name[[2]] = "MVN";          methods_list$args[[2]] = list(blocks = list(block_y = c("Y1","Y2","Y3"), block_x = c("Xcom1")), method = "norm", maxit = 2000, m = 1)
    methods_list$procedure[[3]] = "mice";           methods_list$name[[3]] = "CART";         methods_list$args[[3]] = list(method = "cart", maxit = 100, m = 1)
    methods_list$procedure[[4]] = "mice";           methods_list$name[[4]] = "RF";           methods_list$args[[4]] = list(method = "rf", ntree = 10, maxit = 100, m = 1)
    methods_list$procedure[[5]] = "mice";           methods_list$name[[5]] = "PMM";          methods_list$args[[5]] = list(method = "pmm", maxit = 100, m = 1)
    methods_list$procedure[[6]] = "mice";           methods_list$name[[6]] = "PMM-MT";       methods_list$args[[6]] = list(method = "midastouch", maxit = 100, m = 1)
    
    #### Pre-processing ####
    
    # Identify all conditions to be evaluated across reps = 1:Replications for each method in method_list 
    data_conditions = data.frame(expand.grid(N = N_vec, J_Y = J_Y_vec, J_Xcom = J_Xcom_vec, 
                                             J_Xinc = J_Xinc_vec, MD = MD_vec, class_size = names(pi_list), rho_YX = rho_YX_vec, 
                                             C_modifies_YX = C_modifies_YX_vec, t_rotate = t_rotate_vec))
    
    Kmax = max(simplify2array(lapply(pi_list, FUN = length)))
    
    # Add in the number of classes corresponding to each condition and the marginal probabilities (pi_k)
    pi_df = data.frame(mat.or.vec(nr = length(pi_list), nc = Kmax + 2))+NA
    names(pi_df) = c("class_size","K", paste("pi_", 1:Kmax, sep = ""))
    for(p in 1:nrow(pi_df)){
      pi_df$class_size[p] = names(pi_list)[p]
      pi_df$K[p] = length(pi_list[[p]])
      pi_df[p,seq(3,pi_df$K[p]+2)] = pi_list[[p]]
    }
    data_conditions = merge(x = data_conditions, y = pi_df)
    rm(pi_df)
    
    
    # Ensure that there is at least one data condition
    Ztot = nrow(data_conditions)
    if (Ztot==0){stop("Total number supported of data conditions is zero.")}
    
    
    #Sort and check pctmiss is between 0 and 1
    pctmiss_vec = sort(pctmiss_vec)
    if( (sum(pctmiss_vec<=0) + sum(pctmiss_vec>=1))>0){stop("elements in pctmiss_vec must be between 0 and 1 (exclusive).")}
    
    # Define working directory
    main_wd = paste(dropbox_wd, "/Dissertation/lpa-mi-impute", sep = "")
    
    
    # Check that the dropbox and temp_dir are found on machine
    if(length(temp_dir)==1){
      if (dir.exists(dropbox_wd)*dir.exists(temp_dir)==0){
        stop("Provided dropbox_wd or temp_dir paths not found." )
      }
    } else {
      if (!dir.exists(dropbox_wd)){
        stop("Provided dropbox_wd paths not found." )
      }
      wd_tmp = getwd()
      message("Testing access to each temp_dir:")
      for (p in 1:length(temp_dir)){
        setwd(temp_dir[p])
        message(paste0("  ",temp_dir[p], " successfully accessed"))
        setwd(wd_tmp)
      }
    }
    
    # Create temporary folder
    temp_wd =  paste(temp_dir, temp_fname, sep = "/")
    for (p in 1:length(temp_wd)){
      dir.create(temp_wd[p])
    }
    
    # Create processor-specific temporary folder
    
    if (Processors%%length(temp_wd)!=0){
      stop("Number of working directories is not a multiple of the number of processors. Make sure length(temp_wd)%%Processors == 0.")
    }
    
    if (length(temp_wd)==1){
      setwd(temp_wd)
      for (p in 1:Processors){
        dir.create(paste("p",p, sep = ""))
      }
      temp_wd_p_vec = dir(path = temp_wd, full.names = TRUE)
    } else {
      wt_tmp = getwd()
      temp_wd_p_vec = rep("",Processors)
      n_twd = 1
      for (p in 1:Processors){
        setwd(temp_wd[n_twd])
        dir.create(paste("p",p, sep = ""))
        setwd(paste("p",p, sep = ""))
        
        temp_wd_p_vec[p] = getwd()
        n_twd = ifelse(n_twd==length(temp_wd), 1, n_twd+1)
        setwd(wd_tmp)
      }
      
    } 
    
    
    # Within the processor-specific temporary folder, create subfolders
    for (p in 1:Processors){
      setwd(temp_wd_p_vec[p])
      dir.create("Complete data")
      dir.create("Observed data")
      dir.create("Imputed data")
      dir.create("Stacked data")
      
      for (pm in 1:length(pctmiss_vec)){
        dir.create(paste0(temp_wd_p_vec[p],"/Observed data/pm", pm))
        dir.create(paste0(temp_wd_p_vec[p],"/Imputed data/pm", pm))
        dir.create(paste0(temp_wd_p_vec[p],"/Stacked data/pm", pm))
        
        for (pva in 1:length(methods_list$procedure)){
          dir.create(paste0(temp_wd_p_vec[p],"/Imputed data/pm", pm,"/pva",pva))
          dir.create(paste0(temp_wd_p_vec[p],"/Stacked data/pm", pm,"/pva",pva))
        }
      }
      
    }
    
    
    
    
    # Initiate a list of folder files
    dffolderfiles = NULL
    
    # # writeout the data
    #setwd(main_wd); save.image(file = "image lpa-mi-impute.RData")
    
    # create a list of packages loaded (needed for the "foreach" function)
    temp = strsplit(subset(search(), startsWith(search(),"package")),"package:", fixed = TRUE)
    packages_vec = rep(NA,length(temp))
    for (t in 1:length(temp)){
      packages_vec[t] = temp[[t]][2]
    }
    rm(temp)
    
    
  #### Conduct parallel processing
  #Register processors
  cl<-makePSOCKcluster(Processors)
  registerDoParallel(cl)
    
#p = 1
#rep = 1
#z = 1
  
  # Save environment
  save.image(file = paste0(results_wd,"/environment.RData"))
  
    
  foreach(p = 1:Processors,
          .packages = packages_vec, .options.RNG=42) %dorng%{
        for(z in 1:Ztot){
          pop_params_z = get_FMM_params(z = z, data_conditions = data_conditions)
            for (rep in 1:Replications){
              
              sink(file = paste0(temp_wd_p_vec[p],"/sink p", p, " z", z, " rep", rep, ".txt"))

              dff_complete_rep <- dff_observed_rep <- dff_imputed_rep <- NULL
              
              # Delete any existing files
              print("-----------------Delete any existing Mplus files--------------------------")
              file.remove(list.files(paste0(temp_wd_p_vec[p],"/Observed data"), pattern = ".inp", full.names = TRUE, recursive = TRUE))
              file.remove(list.files(paste0(temp_wd_p_vec[p],"/Observed data"), pattern = ".out", full.names = TRUE, recursive = TRUE))
              file.remove(list.files(paste0(temp_wd_p_vec[p],"/Complete data"), pattern = ".inp", full.names = TRUE, recursive = TRUE))
              file.remove(list.files(paste0(temp_wd_p_vec[p],"/Complete data"), pattern = ".out", full.names = TRUE, recursive = TRUE))
              file.remove(list.files(paste0(temp_wd_p_vec[p],"/Imputed data"), pattern = ".inp", full.names = TRUE, recursive = TRUE))
              file.remove(list.files(paste0(temp_wd_p_vec[p],"/Imputed data"), pattern = ".out", full.names = TRUE, recursive = TRUE))
              
              #### Construct Data Files ####
              print("-----------------Construct .dat Files--------------------------")
              list_complete = lpa.mi.src::get_complete_data(z = z, 
                                                            data_conditions = data_conditions,
                                                            rep = rep, 
                                                            p = p, 
                                                            save_it = TRUE, 
                                                            temp_wd_p_vec = temp_wd_p_vec)
              
              list_observed = lpa.mi.src::get_obs_data(z=z, 
                                                       df =  list_complete$dfcom,
                                                       pctmiss_vec = pctmiss_vec,
                                                       data_conditions = data_conditions,
                                                       save_it = TRUE, 
                                                       temp_wd_p_vec = temp_wd_p_vec,
                                                       p =p , rep = rep)
              
              list_imputed = lpa.mi.src::get_imputed_data(z = z, 
                                                          list_get_obs = list_observed,
                                                          list_get_complete = list_complete,
                                                          methods_list = methods_list,
                                                          data_conditions = data_conditions,
                                                          pctmiss_vec = pctmiss_vec, 
                                                          save_it = TRUE,
                                                          p = p, 
                                                          rep = rep, 
                                                          temp_wd_p_vec = temp_wd_p_vec)
              
              
              # dffilefolders
              dff_complete = list_complete$dffolderfiles
              dff_observed= list_observed$dffolderfiles
              dff_imputed = subset(list_imputed$dffolderfiles, !is.na(m))
              
              # Fit complete model
              print("----------------- Complete Data Results--------------------------")
              out_ftc_complete = fit_til_convergence(p = p , z = z, temp_wd_p_vec = temp_wd_p_vec, 
                                                     dff_target =dff_complete, 
                                                     target_wd = paste0(temp_wd_p_vec[p], "/Complete data"), 
                                                     pop_params_kk = pop_params_z, type_imputation = FALSE, readModels = TRUE,
                                                     savedata = TRUE)
              results_complete = switch_and_get(out_ftc = out_ftc_complete, z = z, data_conditions = data_conditions)$parameters_df %>%
                transform(Data = "Complete",p = p, z = z, rep = rep, pm = NA, pva = NA, m = NA)
              
              # Fit observed data
              print("----------------- Observed Data Results --------------------------")
              results_observed = NULL
              for (pm in 1:length(pctmiss_vec)){
                idx_pm = which(endsWith(as.character(dff_observed$folders),paste0("pm",pm)))
                # Fit the model
                out_ftc_observed = fit_til_convergence(p = p , z = z, temp_wd_p_vec = temp_wd_p_vec, 
                                                       dff_target = dff_observed[idx_pm, ], 
                                                       target_wd = paste0(temp_wd_p_vec[p], "/Observed data/pm",pm),
                                                       pop_params_kk = pop_params_z, type_imputation = FALSE, 
                                                       readModels = TRUE, savedata = TRUE)
                results_observed_pm = switch_and_get(out_ftc = out_ftc_observed, z = z, data_conditions = data_conditions)$parameters_df %>%
                  transform(Data = "Observed",p = p, z = z, rep = rep, pm = pm, pva = NA, m = NA)
                results_observed = rbind(results_observed, results_observed_pm)
              }
              
              
              #Fit the imputed data
              print("----------------- Imputed Data Results--------------------------")
              results_imputed = NULL
              for(pm in 1:length(pctmiss_vec)){ 
                for(pva in 1:length(methods_list$procedure)){
                  M_pva = methods_list$args[[pva]]$m
                  for (m in 1:M_pva){
                    
                    idx_pm_pva_m = which(endsWith(as.character(dff_imputed$folders),paste0("pm",pm,"/pva",pva)) & 
                                           endsWith(as.character(dff_imputed$files), paste0("imp_",m,".dat")))
                    out_ftc_m = fit_til_convergence(p = p , z = z, temp_wd_p_vec = temp_wd_p_vec, 
                                                    dff_target = dff_imputed[idx_pm_pva_m,], 
                                                    target_wd = paste0(temp_wd_p_vec[p], "/Imputed data/pm",pm,"/pva",pva),
                                                    pop_params_kk = pop_params_z, 
                                                    type_imputation = FALSE, m = m, 
                                                    readModels = TRUE, savedata = TRUE)
                    results_imputed_pm_pva_m = switch_and_get(out_ftc = out_ftc_m, z = z, data_conditions = data_conditions)$parameters_df %>%
                      transform(Data = "Imputed",p = p, z = z, rep = rep, pm = pm, pva = pva, m = m)
                    results_imputed = rbind(results_imputed, results_imputed_pm_pva_m)
    
                  } # end for m 
                } # end for pva 
              } # end for pm
              
              
              print("-----------------Write out the results--------------------------")
              
              results_df = rbind(results_complete, results_observed, results_imputed) %>% transform(results_wd = results_wd)
              save(results_df, file = paste0(results_wd,"/p",p," z",z," rep",rep,".RData"))
              
              
              
              
              sink()
            } ## end for rep = 1:Replications
        } #end for z = 1:Ztot 
    
        
    } #for each p = 1,...,P
    
    stopImplicitCluster()
    save.image("environment.RData")
