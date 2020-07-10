  rm(list = ls())
  
  
  library(plyr)
  library(tidyverse)
  library(data.table)
  library(lpa.mi.src)
  library(mice)
  library(dplyr)
  library(doParallel)
  library(foreach)
  library(doRNG)
  require(snow)
  require(doSNOW)
  require(foreach)
  require(pbapply)
  
  computer_name = "MC2"
  Processors = 10
  z_vec = 1:40
  
  # Directories
  dropbox_wd = "D:/Dropbox"
  #dropbox_wd = "C:/Users/marcu/Dropbox"
  results_wd = paste0(dropbox_wd, "/Dissertation/lpa-mi-impute/stage4c-combine-results")
  stage6_wd = paste0(dropbox_wd, "/Dissertation/lpa-mi-impute/stage6-classification-accuracy")
  environment_wd = paste0(dropbox_wd,"/Dissertation/environmental-variables/")
  pingpong_wd = paste0("S:/ping-pong")
  
  system("rm -r H:\\rdata-files")
  system("rm-r H:\\classify-accuracy-files")
  
  Processors = 10
  cl<-makeSOCKcluster(Processors)
  doSNOW::registerDoSNOW(cl)   
  
  
  # Load in the results
  setwd(results_wd)
  load(file ="parameters-combined-results-lpa-mi-impute.RData")
  parameters_combined_df$pva[parameters_combined_df$data_type=="Complete data"] = -1
  parameters_combined_df$pva[parameters_combined_df$data_type=="Observed data"] = 0
  
##Create a parameters cleaned 
#system("rm -r H:\\cleaned-results")
#system("mkdir H:\\cleaned-results")
#system('xcopy "S:\\cleaned-results\\cleaned-parameters-lpa-mi-impute.RData" H:\\cleaned-results /J')
#parameters_cleaned_df = load("H:\\cleaned-results\\cleaned-parameters-lpa-mi-impute.RData")

#rep_x = 467

  for(rep_x in sample(1:500,size = 500, replace = F)){

    if( !(paste0("rep",rep_x,".csv")%in%list.files(path = pingpong_wd)) ){
      
      print(paste0("Replication: ", rep_x))
      tic = proc.time()
      write.csv(x = data.frame(computer = computer_name, total_time = NA), file = paste0(pingpong_wd,"/rep",rep_x,".csv"), row.names = FALSE)    
      
    
            # Make replication directories
            system("rm -r H:\\rdata-files")
            system('mkdir H:\\rdata-files')
            
            system("rm -r H:\\classify-accuracy-files")
            system('mkdir H:\\classify-accuracy-files')
            
            system("rm -r H:\\pooled-trackers")
            system("mkdir H:\\pooled-trackers")
            
            system("rm -r H:\\est-files")
            system("mkdir H:\\est-files")
            
            # Copy over the complete data files
            system(paste0('xcopy "S:\\rdata-files\\list-complete rep',rep_x,' *.RData" H:\\rdata-files'))
            
            # copy over the observed data files
            system(paste0('xcopy "S:\\rdata-files\\list-observed rep',rep_x,' *.RData" H:\\rdata-files'))
            
            # copy over the imputed data files
            system(paste0('xcopy "S:\\rdata-files\\list-imputed rep',rep_x,' *.RData" H:\\rdata-files'))
            
            # copy over the pooled files
            system(paste0('xcopy S:\\pooled-trackers\\pooled-tracker-rep',rep_x,'.RData H:\\pooled-trackers'))
            
            # copy over the est-files
            system(paste0('xcopy S:\\est-files\\est-files-rep',rep_x,'.zip H:\\est-files'))
            
            # unzip
            system(paste0("Bandizip.exe x -y -o:H:\\est-files\\rep",rep_x," H:\\est-files\\est-files-rep",rep_x,".zip"))
            
            
            # load the tracker
            load(paste0("H:/pooled-trackers/pooled-tracker-rep",rep_x,".RData"))
            
            
            
            # create a subpopulation list
            list_subpop<-
              lapply(X = z_vec, FUN = function(zz){
                tmp1 = "complete"; tmp2=".RData";
                load(paste0("H:/rdata-files/list-",tmp1," rep",rep_x," z",zz,tmp2))
                return(list_complete$dfcom %>% select("subpop"))
              })
            
            
            
            out_x = expand.grid(z = z_vec, pva = c(-1:4), pm=1) %>% data.frame() %>% transform(m = NA, data_type = NA) 
            out_x$data_type[out_x$pva==-1] = "Complete data"
            out_x$data_type[out_x$pva==0] = "Observed data"
            out_x$data_type[out_x$pva>0] = "Imputation"
            out_x = out_x %>% transform(kappa1c1=NA, kappa1c2=NA, kappa1c3=NA, kappa2c1=NA, kappa2c2=NA, kappa2c3=NA, kappa3c1=NA, kappa3c2=NA, kappa3c3=NA)
            
            
  
  #x = 1         
            pb <- pbapply::timerProgressBar(max = nrow(out_x), style = 1, width = getOption("width")/4)
            progress <- function(x){setTimerProgressBar(pb, x)}
            opts <- list(progress = progress)
            outlist_x<-
              foreach(x = 1:nrow(out_x),
                      .packages = c("mice","plyr","tidyverse","data.table","dplyr","lpa.mi.src"), 
                      .inorder = TRUE, 
                      .options.snow = opts) %dopar% {
  #for(x in 1:nrow(out_x)){print(x)
                            
                            # Get preliminaries
                            z_x = out_x$z[x]; pva_x = out_x$pva[x]; pm_x = out_x$pm[x]; type_x = out_x$data_type[x]; 
                            
                            # Get the parameters
                            parameters_x = parameters_combined_df %>% 
                              filter(rep==rep_x & z==z_x & pm==pm_x & data_type==type_x & pva==pva_x)
                            
                            if(nrow(parameters_x)>0){
                            
                                if(type_x !="Imputation"){
                                  # Get Q list of observed or complete data
                                    Qlist_x <- parameters_x %>% 
                                      select(paramHeader,param,LatentClass,est) %>% 
                                      Mplus2Qlist()
                                } else {
                                    # Get list of Qlists for imputed data
                                    tracker_x = tracker_df %>% filter(rep==rep_x & z==z_x & pm==pm_x & data_type=="Imputed data" & pva==pva_x & converged == T)
                                    mvec_x = tracker_x$m
                                    mm = 1
                                    Qlist_x = list(NULL)
                                    for (mm in mvec_x){
                                        load(with(tracker_x %>% filter(m == mm),paste0(estwd.x,"/",estfolder.x,"/",estfile)))
                                        Qlist_x[[mm]] = list_estimates$out_Mplus$parameters$unstandardized %>% select(paramHeader,param,LatentClass,est) %>% Mplus2Qlist()
                                    }
                                } #end if(type_x != "Imputation")
                                
                                
                                # Load the data
                                if(type_x=="Complete data"){tmp1 = "complete"; tmp2=".RData"}
                                if(type_x=="Observed data"){tmp1 = "observed"; tmp2=paste0(" pm",pm_x,".RData")}
                                if(type_x=="Imputation"){tmp1 = "imputed"; tmp2 = paste0(" pm",pm_x," pva",pva_x,".RData")}
                                load(paste0("H:/rdata-files/list-",tmp1," rep",rep_x," z",z_x,tmp2))
                                
                                
                                if(type_x!="Imputation"){
                                  return_x = out_x[x, ]
                                  if(type_x=="Complete data"){Y_x = list_complete$dfcom %>% select(starts_with("Y"))}
                                  if(type_x=="Observed data"){Y_x = list_observed$list_obsdf$pm1 %>% select(starts_with("Y"))}
                                  cprob_x <- lpa.mi.src::cprobs(Y_i = Y_x, pi_vec = Qlist_x$pi, mu_mat = Qlist_x$mu, S_array = Qlist_x$S)
                                  modal_x = apply(cprob_x, 1, which.max)
                                  table_x = table(list_subpop[[z_x]]$subpop, modal_x)
                                  return_x$kappa1c1 = table_x[1,1]
                                  return_x$kappa1c2 = table_x[1,2]
                                  return_x$kappa1c3 = table_x[1,3]
                                  return_x$kappa2c1 = table_x[2,1]
                                  return_x$kappa2c2 = table_x[2,2]
                                  return_x$kappa2c3 = table_x[2,3]
                                  return_x$kappa3c1 = table_x[3,1]
                                  return_x$kappa3c2 = table_x[3,2]
                                  return_x$kappa3c3 = table_x[3,3]
                                
                                } else {
                                  tmp_mids = list_imputed$obj_call[[pm_x]][[1]]
                                  return_x<-
                                    lapply(X = mvec_x, 
                                           FUN = function(mm){
                                             Y_x = mice::complete(tmp_mids, action = mm) %>% select(starts_with("Y"));
                                             cprob_x = lpa.mi.src::cprobs(Y_i = Y_x, pi_vec = Qlist_x[[mm]]$pi, mu_mat = Qlist_x[[mm]]$mu, S_array = Qlist_x[[mm]]$S) %>% data.frame() %>% transform(id = 1:nrow(Y_x), m = mm)
                                             modal_x = apply(cprob_x %>% select(starts_with("x")), 1, which.max)
                                             table_x = table(list_subpop[[z_x]]$subpop, modal_x)
                                             return(data.frame(z=z_x, pva=pva_x, pm=pm_x, m=mm, data_type="Imputation",
                                                               kappa1c1 = table_x[1,1], kappa1c2=table_x[1,2], kappa1c3=table_x[1,3], 
                                                               kappa2c1=table_x[2,1], kappa2c2=table_x[2,2], kappa2c3=table_x[2,3], 
                                                               kappa3c1=table_x[3,1], kappa3c2=table_x[3,2], kappa3c3=table_x[3,3]))
                                           } ) %>% data.table::rbindlist() %>% arrange(z,pva,m)
                                }
                                
                                
                                return(return_x)
                
                } #if(nrow(parameters_x)>0)
              }#end for(x = )
            
            out_x<-rbindlist(outlist_x) %>% data.frame()
            
            out_x$pva[out_x$data_type!="Imputation"] = NA
            save(out_x, file = paste0("H:/classify-accuracy-files/classify-accuracy-rep",rep_x,".RData"))
            system(paste0('xcopy H:\\classify-accuracy-files\\classify-accuracy-rep',rep_x,'.RData S:\\classify-accuracy-files /J /Y'))
  
            toc = proc.time()-tic; toc = round(toc[[3]],0)
            write.csv(x = data.frame(computer = computer_name, total_time = toc), file = paste0(pingpong_wd,"/rep",rep_x,".csv"), row.names = FALSE)    
            
            # Clean up
            system("rm -r H:\\rdata-files")
            system("rm -r H:\\classify-accuracy-files")
            system("rm -r H:\\pooled-trackers")
            system("rm -r H:\\est-files")
  
       }#end if(ping)
   }# end for rep=

   stopCluster(cl)

