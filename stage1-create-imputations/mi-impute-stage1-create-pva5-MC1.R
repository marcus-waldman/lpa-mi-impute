#*****NOTE THAT THIS CODE ASSUMES THAT STAGE0 LISTS ARE EXTRACTED AND READY TO GO*************

    #  .rs.restartR()
    
    rm(list = ls())
    
    library(plyr)
    library(tidyverse)
    library(lpa.mi.src)
    library(doParallel)
    library(foreach)
    library(doRNG)
    require(snow)
    require(doSNOW)
    require(foreach)
    require(pbapply)
    library(stringr)
    
    computer_name = "MC1"
    pingpong_wd = "S:/ping-pong"
    
    #
    #
    #
    # system('mkdir H:\\Rdata-mi-impute-stage0')
    # system('xcopy "S:\\Rdata-mi-impute-stage0 Jan 06 2019 13 29.zip" H:\\Rdata-mi-impute-stage0 /Y /J')
    # system('Bandizip x -y "H:\\Rdata-mi-impute-stage0\\Rdata-mi-impute-stage0 Jan 06 2019 13 29.zip"')
    # system('del "H:\\Rdata-mi-impute-stage0\\Rdata-mi-impute-stage0 Jan 06 2019 13 29.zip"')
    # 
    # system('mkdir H:\\EMs-temp')
    # system('mkdir E:\\dat-files-pva5')
    
    
    ######################################################################################################################
    #      Load in and select the data generating conditions                                                             #
    ######################################################################################################################
    # Gaby
    load("S:/environment-mi-impute-stage0 Jan 06 2019 13 29.RData")    
    
    
    
    ######################################################################################################################
    #      (Re)Define some simulation variables and such                                                                 #
    ######################################################################################################################
    Processors = 10
    this_file = "mi-impute-stage1-create-pva5-QMR.R"
    
    # Acer-laptop
      #dropbox_wd =  "C:/Users/marcu/Dropbox"  #  Dropbox folder on local machine
      #datfiles_dir = "C:/Users/marcu/Documents/lpa-mi-create-stage0" #Stores Mplus-read dat files
      #Rdata_dir = "C:/Users/marcu/Documents/lpa-mi-create-stage0"
    # Queen Mary's Revenge
      #dropbox_wd =  "D:/Dropbox"  #  Dropbox folder on local machine
      #datfiles_dir = "C:/Users/marcu/Documents" #Stores Mplus-read dat files
      #Rdata_dir = "C:/Users/marcu/Documents"  
    # Xi-GSU
      #dropbox_wd = "C:/Users/mwaldman1/Dropbox"
      #datfiles_dir =  #Stores Mplus-read dat files
    # Big Bertha
      #dropbox_wd =  "D:/Dropbox/Dropbox"  #  Dropbox folder on local machine
      #datfiles_dir = #Stores Mplus-read dat files
      #Rdata_dir = 
    # Gaby/MC/QMR
      dropbox_wd = "D:/Dropbox"
      datfiles_dir = "E:" #Stores Mplus-read dat files
      Rdata_dir = "E:"
      temp_dir = "H:"
      tempdir_EMs = "H:/EMs-temp"
      stage0_dir = "H:/Rdata-mi-impute-stage0"
      

    # Define the imputation alogrithms to compare
      M = 100
      methods_list = list(procedure = list(NULL), name = list(NULL), args = list(NULL))
      pva_vec = c(5)
      methods_list$procedure[[1]] = "EMs_LPA"; methods_list$name[[1]] = "EM-Sampling";    methods_list$args[[1]] = list(M = M, tempdir_mplus = tempdir_EMs, starts0 = 20, bayes_boot = TRUE, boot_2x = FALSE)
      methods_list$pva_vec = pva_vec
      
      
    ######################################################################################################################
    #      Select on desired data conditions                                                                             #
    ######################################################################################################################
    data_conditions = transform(data_conditions, z = 1:nrow(data_conditions))
    zvec =1:nrow(data_conditions)
      

    ######################################################################################################################
    #      Create folders and subfolders in temporary directory for read-ins and write-outs at each replication          #
    ######################################################################################################################
     timestamp = format(Sys.time(), "%b %d %Y %H %M")
     # Create folder name for temporary folder (saved in results and in temp files)
      temp_fname = paste("mi-impute-stage1-create", timestamp , sep = " ")
      # Create temporary folder
      temp_wd =  paste(temp_dir, temp_fname, sep = "/")
      dir.create(temp_wd)
    


    ######################################################################################################################
    #    Create file that stores results at each replication eventually to be zipped and copied to dropbox results       #
    ######################################################################################################################
    lists_wd = paste0(Rdata_dir,"/Rdata-",temp_fname)
    dir.create(lists_wd)
    # Save environment
    save.image(file = paste0(lists_wd,"/environment-", temp_fname,".RData"))
    # Copy all R files in main directory
    file.copy(from = paste0(main_wd,"/",this_file), to = paste0(lists_wd,"/","executed-R-code ", temp_fname, ".R"))
    
    
    
    
    ######################################################################################################################
    #      Initiate other variables and write out the environment                                                        #
    ######################################################################################################################
    # Set up parallel processing 
    cl<-makeSOCKcluster(Processors)
    registerDoSNOW(cl)

#    rep = 1
#    z = zvec[1]
    


    for(rep in 1:Replications){
      if( !(paste0("rep",rep,".csv")%in%list.files(path = pingpong_wd)) ){
        
        tic = proc.time()
        write.csv(x = data.frame(computer = computer_name, total_time = NA), file = paste0(pingpong_wd,"/rep",rep,".csv"), row.names = FALSE)    
        
            # Create replication-specific temporary folder
            datfiles_wd_rep_vec = paste0(temp_wd,"/rep",1:500, sep = "")
            dir.create(datfiles_wd_rep_vec[rep])
            #}
      
            # Within the processor-specific temporary folder, create subfolders
            #for (rep in 1:Replications){
            for(z in zvec){
              for (ii in 1:length(pctmiss_vec)){
                for (jj in 1:length(methods_list$procedure)){
                  dir.create(paste0(datfiles_wd_rep_vec[rep],"/z",z,"/Imputed data/pm", pm_vec[ii],"/pva",methods_list$pva_vec[[jj]]), recursive = TRUE)
                } #for ii
              } # for jj
            } #for z
            #} #for rep
            
            
            
            
            for(z in zvec){
              load(paste0(stage0_dir, "/list-complete rep",rep, " z", z,".RData"))
              for(ii in 1:length(pctmiss_vec)){
                load(paste0(stage0_dir, "/list-observed rep",rep, " z", z," pm", pm_vec[ii],".RData"))
                
                
                list_imputed = get_imputed_data(z = z, list_get_obs = list_observed, list_get_complete = list_complete, methods_list = methods_list, 
                                                data_conditions = data_conditions, pctmiss_vec = pctmiss_vec, save_it = TRUE, rep = rep, cl_obj = cl, 
                                                temp_wd_rep_vec = datfiles_wd_rep_vec)
                
                
                fname_imputed = paste0(lists_wd, "/list-imputed rep", rep, " z", z, " pm", pm_vec[ii], " pva", pva_vec[jj],".RData")
                #save(list_imputed, file = fname_imputed)
                save(list_imputed, file = paste0("S:/Rdata-mi-impute-pva5/list-imputed rep", rep, " z", z, " pm", pm_vec[ii], " pva", pva_vec[jj],".RData"))
                print(paste0(">>list_imputed: ",fname_imputed," saved"))
                print(" ")
                
                
                rm(list_observed)
                rm(list_imputed)
      
              } #end for ii
              rm(list_complete)
            } ## end for z = 1:Ztot
      #      setWinProgressBar(pb, rep, title=paste( round(rep/Replications*100, 0),"% done"))
          
          print("Zipping dat files and saving...")
          setwd(datfiles_wd_rep_vec[rep]) 
          setwd("..")
          system(paste0('Bandizip c -y -l:2 E:\\dat-files-pva5\\dat-files-pva5-rep',rep,'.zip rep', rep))
          system(paste0('xcopy  E:\\dat-files-pva5\\dat-files-pva5-rep',rep,'.zip S:\\dat-files-pva5 /J /Y'))
          
          # Clean up
          
    
        # write out estimation time
        toc = proc.time()-tic; toc = round(toc[[3]],0)
        write.csv(x = data.frame(computer = computer_name, total_time = toc), file = paste0(pingpong_wd,"/rep",rep,".csv"), row.names = FALSE)    
        print(" ")
        print(" ")
          
        
      } #End if
          
    } #for rep in 1:Replications


stopCluster(cl)
    

    
