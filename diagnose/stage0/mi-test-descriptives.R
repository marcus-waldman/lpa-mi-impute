#  .rs.restartR()
    
    rm(list = ls())
    
    library(plyr)
    library(tidyverse)
    library(lpa.mi.src)
    library(semTools)

    
    Replications = 500
    stage0_zip = "Rdata-mi-impute-stage0 Jan 02 2019 18 24.zip"

    # Acer-laptop
      dropbox_wd =  "C:/Users/marcu/Dropbox"  #  Dropbox folder on local machine
    # Big Bertha
      #dropbox_wd =  "D:/Dropbox/Dropbox"  #  Dropbox folder on local machine

    # Fully-crossed simulation conditions
      pctmiss_vec = c(0.5); pm_vec = c(1)   #  Percent of observations with *at least* one missing value (indicator or missing data correlate) -> if other specification needed one must modify list_inds.miss generation in get_obs_data.R
      N_vec = c(300, 1200)                  #  Sample sizes within each replication (small - 300, large - 1200)
      J_Y_vec = c(4)                        #  Number of latent class indicators (K = 4)
      J_Xcom_vec = c(1)                     #  Number of complete data missing data correlates
      J_Xinc_vec = c(0)                     #  Number of missing data correlates which themselves contain missing data
      MD_vec = c(2.87,3.7)                  #  Mahalanobis disance separating classes (low - 2.87, high = 3.70)
      pi_list = list(                       #  Mixing and (by proxy) the number of classes, K (maximum of 3)
        A = rep(1/3,3), 
        B = c(0.45,0.45,0.1) 
      )
      rho_YX_vec = c(0.4, 0)              #  Strength of missing data correlates, correlation
      dX_vec = c(TRUE, FALSE)             #  Between class differences in missing data correlates, Cohen's d
      C_modifies_YX_vec = c(FALSE, TRUE)  #  Whether or not the Y-X relationships are different between latent classes 
      t_rotate_vec = c(0)                 #  Angle for rotation matrix for LPA indicators 
      M = 5 #  Number of imputations

    # Define the imputation alogrithms to compare
      methods_list = list(procedure = list(NULL), name = list(NULL), args = list(NULL))
      pva_vec = c(1,7)
      methods_list$procedure[[1]] = "amelia";         methods_list$name[[1]] = "mg-Amelia";    methods_list$args[[1]] = list(m = M, p2s = as.numeric(VerboseImpute), tolerance = 1E-4, cs = "subpop", intercs = TRUE, empri = 0.01)
      #methods_list$procedure[[2]] = "amelia";         methods_list$name[[2]] = "Amelia";       methods_list$args[[2]] = list(m = M, p2s = as.numeric(VerboseImpute), tolerance = 1E-4, empri = 0.01)
      #methods_list$procedure[[3]] = "mice";           methods_list$name[[3]] = "PMM";          methods_list$args[[3]] = list(method = "pmm", maxit = 40, m = M, printFlag = VerboseImpute)
      #methods_list$procedure[[4]] = "mice";           methods_list$name[[4]] = "CART";         methods_list$args[[4]] = list(method = "cart", maxit = 40, m = M, printFlag = VerboseImpute)
      #methods_list$procedure[[5]] = "mice";           methods_list$name[[5]] = "RF";           methods_list$args[[5]] = list(method = "rf", ntree = 10, maxit = 20, m = M, printFlag = VerboseImpute)
      #methods_list$procedure[[6]] = "mice";           methods_list$name[[6]] = "MVN";          methods_list$args[[6]] = list(blocks = list(block_y = c("Y1","Y2","Y3"), block_x = c("Xcom1")), method = "norm", maxit = 5E1, m = M, printFlag = VerboseImpute)
      methods_list$procedure[[2]] = "stratamelia";         methods_list$name[[2]] = "strata-Amelia";    methods_list$args[[2]] = list(m = M, p2s = as.numeric(VerboseImpute), tolerance = 1E-4)
      methods_list$pva_vec = pva_vec
      
      
      
      
    ######################################################################################################################
    #      Pre-process the data generating conditions                                                                    #
    ######################################################################################################################
    # Identify all conditions to be evaluated across reps = 1:Replications for each method in method_list 
      data_conditions = data.frame(expand.grid(N = N_vec, J_Y = J_Y_vec, J_Xcom = J_Xcom_vec, 
                                               J_Xinc = J_Xinc_vec, MD = MD_vec, class_size = names(pi_list), 
                                               rho_YX = rho_YX_vec, dX = dX_vec,
                                               C_modifies_YX = C_modifies_YX_vec, t_rotate = t_rotate_vec))
    # Remove the AV as non-confounder conditions
      z_remove = with(data_conditions, which(rho_YX==0 & dX == FALSE | rho_YX==0 & C_modifies_YX_vec==TRUE))
      if(length(z_remove)>0){data_conditions = data_conditions[-z_remove, ]}
    # Add in the number of classes corresponding to each condition and the marginal probabilities (pi_k)
      Kmax = max(simplify2array(lapply(pi_list, FUN = length)))
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
    # Initiate a list of folder files
      dffolderfiles = NULL
      
      
    ######################################################################################################################
    #      Create folders and subfolders in temporary directory for read-ins and write-outs at each replication          #
    ######################################################################################################################
    # Create folder name for temporary folder (saved in results and in temp files)
      temp_fname = paste("mi-impute-test-mgAmelia", format(Sys.time(), "%b %d %Y %H %M"), sep = " ")
    # Create temporary folder
      temp_wd =  paste(temp_dir, temp_fname, sep = "/")
      dir.create(temp_wd)
    # Unzip the file 
      stage0_dir = paste0(dropbox_wd,"/Dissertation/data-manipulation/lists-stage0")
      dir.create(paste0(temp_wd,"/stage0"))
    # Create replication-specific temporary folder
      datfiles_wd_rep_vec = NULL
      for (rep in 1:Replications){
        temp_wd_rep = paste0(temp_wd,"/rep",rep, sep = "")
        dir.create(temp_wd_rep)
        datfiles_wd_rep_vec = c(datfiles_wd_rep_vec, temp_wd_rep)
      }
      rm(temp_wd_rep)

    # Within the processor-specific temporary folder, create subfolders
      for (rep in 1:Processors){
        setwd(datfiles_wd_rep_vec[rep])
        dir.create("Imputed data")
        for (ii in 1:length(pctmiss_vec)){
          dir.create(paste0(datfiles_wd_rep_vec[rep],"/Imputed data/pm", pm_vec[ii]))
          for (jj in 1:length(methods_list$procedure)){
            dir.create(paste0(datfiles_wd_rep_vec[rep],"/Imputed data/pm", pm_vec[ii],"/pva",methods_list$pva_vec[[jj]]))
          } #for pva
        } # for pm 
      } #for rep 
      
      
      
      ######################################################################################################################
      #    Create file that stores results at each replication eventually to be zipped and copied to dropbox results       #
      ######################################################################################################################
      lists_wd = paste0(Rdata_dir,"/Rdata-",temp_fname)
      dir.create(lists_wd)
      # Save environment
      save.image(file = paste0(lists_wd,"/environment.RData"))
      # Copy all R files in main directory
      file.copy(from = paste0(main_wd,"/",this_file), to = paste0(lists_wd,"/","executed-R-code ", temp_fname, ".R"))
      
       
   ######################################################################################################################
   #      Initiate other variables and write out the environment                                                        #
   ######################################################################################################################
    t1 = proc.time()      
    foreach(rep = 1:Replications, .packages = packages_vec, .options.RNG=RndSeed, .errorhandling = "pass", .inorder = TRUE) %dorng%{
        sink(file = paste0(lists_wd,"/sink rep", rep, ".txt"))
    
            for(z in 1:Ztot){
                 print(paste0("--------------- rep = ", rep, ",  z = ", z, "  ---------------"))
              
                  pop_params_z = get_FMM_params(z = z, data_conditions = data_conditions)
                  
                  
                  # Unzip and load complete and observed data
                  file_unzip = paste0("list-stage0 rep", rep, " z",z,".RData")
                  dir_unzip = paste0(temp_wd,"/stage0")
                  unzip(zipfile = paste(stage0_dir, stage0_zip, sep = "/"), 
                        files = file_unzip, 
                        exdir = dir_unzip)
                  load(paste0(dir_unzip,"/",file_unzip))
                  
                  # Generate imputations
                  list_imputed = lpa.mi.src::get_imputed_data(z = z, 
                                                              list_get_obs = list_observed,
                                                              list_get_complete = list_complete,
                                                              methods_list = methods_list,
                                                              data_conditions = data_conditions,
                                                              pctmiss_vec = pctmiss_vec,
                                                              save_it = FALSE,
                                                              rep = rep, 
                                                              temp_wd_rep_vec = datfiles_wd_rep_vec)
    #              print("  >> list_imputed constructed")
                  
                  # Fit the a model
                  obj_sem = sem(model = form_lavaan, 
                                         data = list_complete$dfcom, 
                                         group = "subpop", 
                                         estimator = "ml", information = "observed")
                  
                  obj_smi_pva1 = sem.mi(model = form_lavaan, 
                                data = mids2datlist(list_imputed$obj_call$pm1$pva1),
                                group = "subpop",
                                estimator = "ml", information = "observed")
                  
                  obj_smi_pva7 = sem.mi(model = form_lavaan, 
                                       data = mids2datlist(list_imputed$obj_call$pm1$pva7),
                                       group = "subpop",
                                       estimator = "ml", information = "observed")
                  
                  
                  # Save the model results
                  vars_keep = c("lhs","op","rhs","group","est","se")
                  df_complete = obj_sem %>% parameterEstimates() %>% select(vars_keep) %>% transform(rep=rep,z=z,pva=NA,type="Complete")
                  df_mga = obj_smi_pva1 %>% summary() %>% data.frame() %>% select(vars_keep) %>% transform(rep=rep,z=z,pva=1,type="Imputed")
                  df_sta = obj_smi_pva7 %>% summary() %>% data.frame() %>% select(vars_keep) %>% transform(rep=rep,z=z,pva=7,type="Imputed")
                  
                  out_df = rbind(df_complete, df_mga, df_sta)  
                  save(out_df, file = paste0(lists_wd,"/result-test-mgAmelia rep",rep, " z",z, " pva",1, ".rds"))
                  
                  rm(list_complete)
                  rm(list_observed)
                  rm(list_imputed)
                  
    
               } ## end for z = 1:Ztot
        sink()
        } #foreach rep = 1,...,Replications
    stopImplicitCluster()
    
    proc.time()-t1    
        

              
