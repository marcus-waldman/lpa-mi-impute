#  .rs.restartR()
    
    rm(list = ls())
    
    library(plyr)
    library(tidyverse)
    library(lpa.mi.src)
    library(doParallel)
    library(foreach)
    library(doRNG)

    
    Processors = 55
    Replications = 10
    RndSeed = 42
    VerboseImpute = FALSE
    this_file = "mi-compare-existing.R"
    
    # Acer-laptop
      #dropbox_wd =  "C:/Users/marcu/Dropbox"  #  Dropbox folder on local machine
      #temp_dir = "C:/Users/marcu/Documents/lpa-mi-impute" #  Temporary directory for local maching (will automaticall add new temporary folder)
    # Queen Mary's Revenge
      #dropbox_wd =  "D:/Dropbox"  #  Dropbox folder on local machine
      #temp_dir = "C:/Users/marcu/Documents/lpa-mi-impute" #  Temporary directory for local maching (will automaticall add new temporary folder)        
    # Xi-GSU
      #dropbox_wd = "C:/Users/mwaldman1/Dropbox"
      #temp_dir = "C:/Users/mwaldman1/Documents/lpa-mi-impute"
    # Big Bertha
      #dropbox_wd =  "D:/Dropbox/Dropbox"  #  Dropbox folder on local machine
      #temp_dir = "D:/lpa-mi-impute" #  Temporary directory for local maching (will automaticall add new temporary folder)        
    # Gaby
      dropbox_wd = "Z:/Dropbox"
      temp_dir = "C:/Users/Gaby/Documents/lpa-mi-src" 
      zip_dir = "C:/Users/Gaby/Documents/zip-temp" 


    # Fully-crossed simulation conditions
      pctmiss_vec = c(0.5)          #  Percent of observations with *at least* one missing value (indicator or missing data correlate) -> if other specification needed one must modify list_inds.miss generation in get_obs_data.R
      N_vec = c(300, 1200)          #  Sample sizes within each replication (small - 300, large - 1200)
      J_Y_vec = c(4)                #  Number of latent class indicators (K = 4)
      J_Xcom_vec = c(1)             #  Number of complete data missing data correlates
      J_Xinc_vec = c(0)             #  Number of missing data correlates which themselves contain missing data
      MD_vec = c(2.87,3.7)          #  Mahalanobis disance separating classes (low - 2.87, high = 3.70)
      pi_list = list(               #  Mixing and (by proxy) the number of classes, K (maximum of 3)
        A = rep(1/3,3), 
        B = c(0.45,0.45,0.1) 
      )
      rho_YX_vec = c(0.4, 0)          #  Strength of missing data correlates, correlation
      dX_vec = c(TRUE, FALSE)         #  Between class differences in missing data correlates, Cohen's d
      C_modifies_YX_vec = c(FALSE)    #  Whether or not the Y-X relationships are different between latent classes 
      t_rotate_vec = c(0)             #  Angle for rotation matrix for LPA indicators 
      M = 50 #  Number of imputations

    # Define the imputation alogrithms to compare
      methods_list = list(procedure = list(NULL), name = list(NULL), args = list(NULL))
      methods_list$procedure[[1]] = "stratamelia";    methods_list$name[[1]] = "mg-Amelia";    methods_list$args[[1]] = list(m = M, p2s = as.numeric(VerboseImpute), tolerance = 1E-4)
      methods_list$procedure[[2]] = "amelia";         methods_list$name[[2]] = "Amelia";       methods_list$args[[2]] = list(m = M, p2s = as.numeric(VerboseImpute), tolerance = 1E-4)
      methods_list$procedure[[3]] = "mice";           methods_list$name[[3]] = "CART";         methods_list$args[[3]] = list(method = "cart", maxit = 50, m = M, printFlag = VerboseImpute)
      #methods_list$procedure[[4]] = "mice";           methods_list$name[[4]] = "PMM";          methods_list$args[[4]] = list(method = "pmm", maxit = 50, m = M, printFlag = VerboseImpute)
      #methods_list$procedure[[5]] = "mice";           methods_list$name[[5]] = "RF";           methods_list$args[[5]] = list(method = "rf", ntree = 10, maxit = 20, m = M, printFlag = VerboseImpute)
      #methods_list$procedure[[6]] = "mice";           methods_list$name[[6]] = "PMM-MT";       methods_list$args[[6]] = list(method = "midastouch", maxit = 20, m = M,  printFlag = VerboseImpute)
      #methods_list$procedure[[2]] = "mice";           methods_list$name[[2]] = "MVN";          methods_list$args[[2]] = list(blocks = list(block_y = c("Y1","Y2","Y3"), block_x = c("Xcom1")), method = "norm", maxit = 1E3, m = M, printFlag = VerboseImpute)
      
      
      
    ######################################################################################################################
    #      Pre-process the data generating conditions                                                                    #
    ######################################################################################################################
    # Identify all conditions to be evaluated across reps = 1:Replications for each method in method_list 
      data_conditions = data.frame(expand.grid(N = N_vec, J_Y = J_Y_vec, J_Xcom = J_Xcom_vec, 
                                               J_Xinc = J_Xinc_vec, MD = MD_vec, class_size = names(pi_list), 
                                               rho_YX = rho_YX_vec, dX = dX_vec,
                                               C_modifies_YX = C_modifies_YX_vec, t_rotate = t_rotate_vec))
    # Remove the AV as non-confounder conditions
      z_remove = with(data_conditions, which(rho_YX==0 & dX == FALSE))
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
      temp_fname = paste("sim-lpa-mi-impute", format(Sys.time(), "%b %d %Y %H %M"), sep = " ")
    # Create temporary folder
      temp_wd =  paste(temp_dir, temp_fname, sep = "/")
      for (p in 1:length(temp_wd)){
        dir.create(temp_wd[p])
      }
    # Create processor-specific temporary folder
#      if (Processors%%length(temp_wd)!=0){
#        stop("Number of working directories is not a multiple of the number of processors. Make sure length(temp_wd)%%Processors == 0.")
#      }
      
      if (length(temp_wd)==1){
        setwd(temp_wd)
        for (p in 1:Processors){
          dir.create(paste("p",p, sep = ""))
        }
        temp_wd_p_vec = dir(path = temp_wd, full.names = TRUE)
      } else {
        wd_tmp = getwd()
        temp_wd_p_vec = rep("",Processors)
        n_twd = 1
        for (p in 1:Processors){
          setwd(temp_wd[n_twd])
          dir.create(paste("p",p, sep = ""))
          setwd(paste("p",p, sep = ""))
          
          temp_wd_p_vec[p] = getwd()
          n_twd = ifelse(n_twd==length(temp_wd), 1, n_twd+1)
          setwd(wd_tmp)
        } # end for p
      } # end if (length(temp_wd))
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
      

   ######################################################################################################################
   #           Initate for parallel processing                                                                          #     
   ######################################################################################################################
     # create a list of packages loaded (needed for the "foreach" function)
     temp = strsplit(subset(search(), startsWith(search(),"package")),"package:", fixed = TRUE)
     packages_vec = rep(NA,length(temp))
     for (t in 1:length(temp)){
       packages_vec[t] = temp[[t]][2]
     }
     rm(temp)
     #Register processors
     cl<-makePSOCKcluster(Processors)
     registerDoParallel(cl)

  

   ######################################################################################################################
   #    Create file that stores results at each replication eventually to be zipped and copied to dropbox results       #
   ######################################################################################################################
     tozip_wd = paste0(zip_dir,"/zip-",temp_fname)
     dir.create(tozip_wd)
    # Save environment
       save.image(file = paste0(tozip_wd,"/environment.RData"))
    # Copy all R files in main directory
       file.copy(from = paste0(main_wd,"/",this_file), to = paste0(tozip_wd,"/","executed-R-code ", temp_fname, ".R"))
       
      
       
   ######################################################################################################################
   #      Initiate other variables and write out the environment                                                        #
   ######################################################################################################################
      
      
#p = 1
#rep = 1
#z = 1
  foreach(p = 1:Processors, .packages = packages_vec, .options.RNG=RndSeed, .errorhandling = "pass") %dorng%{
      for (rep in 1:Replications){ 
        for(z in 1:Ztot){
              pop_params_z = get_FMM_params(z = z, data_conditions = data_conditions)
              
#              sink(file = paste0(temp_wd_p_vec[p],"/sink p", p, " z", z, " rep", rep, ".txt"))

              dff_complete_rep <- dff_observed_rep <- dff_imputed_rep <- NULL
              
              # Delete any existing files
              print("-----------------Delete any existing Mplus files--------------------------")
              file.remove(list.files(paste0(temp_wd_p_vec[p],"/Observed data"), pattern = ".inp", full.names = TRUE, recursive = TRUE))
              file.remove(list.files(paste0(temp_wd_p_vec[p],"/Observed data"), pattern = ".out", full.names = TRUE, recursive = TRUE))
              file.remove(list.files(paste0(temp_wd_p_vec[p],"/Observed data"), pattern = ".dat", full.names = TRUE, recursive = TRUE))
              file.remove(list.files(paste0(temp_wd_p_vec[p],"/Complete data"), pattern = ".inp", full.names = TRUE, recursive = TRUE))
              file.remove(list.files(paste0(temp_wd_p_vec[p],"/Complete data"), pattern = ".out", full.names = TRUE, recursive = TRUE))
              file.remove(list.files(paste0(temp_wd_p_vec[p],"/Complete data"), pattern = ".dat", full.names = TRUE, recursive = TRUE))
              file.remove(list.files(paste0(temp_wd_p_vec[p],"/Imputed data"), pattern = ".inp", full.names = TRUE, recursive = TRUE))
              file.remove(list.files(paste0(temp_wd_p_vec[p],"/Imputed data"), pattern = ".out", full.names = TRUE, recursive = TRUE))
              file.remove(list.files(paste0(temp_wd_p_vec[p],"/Imputed data"), pattern = ".dat", full.names = TRUE, recursive = TRUE))
              
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
              
              # Save the descriptives
              #save(list_complete, file = paste0(tozip_wd,"/list_complete p",p," z", z, " rep", rep,".RData"))
              #save(list_observed, file = paste0(tozip_wd,"/list_observed p",p," z", z, " rep", rep,".RData"))
              
              
              # dffilefolders
              dff_complete = list_complete$dffolderfiles
              dff_observed= list_observed$dffolderfiles
              dff_imputed = subset(list_imputed$dffolderfiles, !is.na(m))
             
              # Fit and save the complete model
              print("----------------- Complete Data Results--------------------------")
              out_ftc_complete = fit_til_convergence(p = p , z = z, temp_wd_p_vec = temp_wd_p_vec, 
                                                     dff_target =dff_complete, 
                                                     target_wd = paste0(temp_wd_p_vec[p], "/Complete data"), 
                                                     pop_params_kk = pop_params_z, type_imputation = FALSE, readModels = TRUE,
                                                     savedata = TRUE, output_txt = "tech3")
              switched_complete = resolve_label_switch(out_ftc = out_ftc_complete)
              outlist_complete = to_outlist(out_ftc = out_ftc_complete, out_switched = switched_complete, Data = "Complete", 
                                            p = p, z=z, rep=rep)
              save(outlist_complete, file = paste0(tozip_wd,"/results-list-complete p",p," z", z, " rep", rep,".RData"))
              
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
                                                       readModels = TRUE, savedata = TRUE, output_txt = "tech3")
                switched_observed_pm = resolve_label_switch(out_ftc = out_ftc_observed)
                outlist_observed_pm = to_outlist(out_ftc = out_ftc_observed, out_switched = switched_observed_pm, Data = "Observed", 
                                              p = p, z=z, rep=rep, pm=pm)
               save(outlist_observed_pm, file = paste0(tozip_wd,"/results-list-observed p",p," z", z, " rep", rep," pm",pm,".RData"))
              } # end pm in 1:length(pctmiss_vec)

              #Fit the imputed data
              print("----------------- Imputed Data Results--------------------------")
              for(pm in 1:length(pctmiss_vec)){
                for(pva in 1:length(methods_list$procedure)){

                  out_pool = fit_and_do_rubinrules(dff_imputed = dff_imputed, methods_list = methods_list, pop_params_z = pop_params_z,
                                          data_conditions = data_conditions, temp_wd_p_vec = temp_wd_p_vec,
                                          p = p, z = z, pm = pm, pva = pva)
                  out_pool$Data = "Imputed"; out_pool$p=p; out_pool$z=z; out_pool$rep=rep; out_pool$pm=pm; out_pool$pva=pva
                 save(out_pool, file = paste0(tozip_wd,"/results-list-imputed p",p," z", z, " rep", rep," pm",pm," pva",pva,".RData"))                  
                } # end for pva
              } # end for pm

#             sink()
           } ## end for z = 1:Ztot
       } #end for rep = 1:Replications 
    } #for each p = 1,...,P
    
    stopImplicitCluster()
    
    
    
    ######################################################################################################################
    #    Create file that stores results at each replication eventually to be zipped and copied to dropbox results       #
    ######################################################################################################################
    setwd(tozip_wd)
    files_tozip = list.files(path = tozip_wd, full.names = FALSE, pattern = ".RData")
    simulation_list = sapply(files_tozip,function(x){l=list(mget(load(x))); names(l) =x; return(l)})
    save(simulation_list, file = paste0(dropbox_wd,"/Dissertation/results/",temp_fname,".Rdata"))
