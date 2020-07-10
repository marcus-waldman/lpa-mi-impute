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
  
  computer_name = "MC4"
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
  
  #Load the data conditions
  setwd(environment_wd)
  load(file = "data-conditions.RData")
  
  # Load in the results
  setwd(results_wd)
  load(file ="parameters-combined-results-lpa-mi-impute.RData")
  parameters_combined_df$pva[parameters_combined_df$data_type=="Complete data"] = -1
  parameters_combined_df$pva[parameters_combined_df$data_type=="Observed data"] = 0


  
  
#Create a parameters cleaned 
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
            
            # Copy over the complete data files
            system(paste0('xcopy "S:\\rdata-files\\list-complete rep',rep_x,' *.RData" H:\\rdata-files'))
            
            # copy over the observed data files
            system(paste0('xcopy "S:\\rdata-files\\list-observed rep',rep_x,' *.RData" H:\\rdata-files'))
            
            # copy over the imputed data files
            system(paste0('xcopy "S:\\rdata-files\\list-imputed rep',rep_x,' *.RData" H:\\rdata-files'))
            
            #system("rm -r H:\\est-files")
            #system("mkdir H:\\est-files")
            #
            ## copy over the pooled files
            #system(paste0('xcopy S:\\pooled-trackers\\pooled-tracker-rep',rep_x,'.RData H:\\pooled-trackers'))
            ## copy over the est-files
            #system(paste0('xcopy S:\\est-files\\est-files-rep',rep_x,'.zip H:\\est-files'))
            #
            # unzip
            #system(paste0("Bandizip.exe x -y -o:H:\\est-files\\rep",rep_x," H:\\est-files\\est-files-rep",rep_x,".zip"))
            #
            ## load the tracker
            #load(paste0("H:/pooled-trackers/pooled-tracker-rep",rep_x,".RData"))
            
            
            
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
            out_x$summary_type = NA
            out_x = out_x %>% transform(kl_cprob=NA, bias_cprob1=NA, bias_cprob2=NA, bias_cprob3=NA, 
                                                     rmse_cprob1=NA, rmse_cprob2=NA, rmse_cprob3=NA, 
                                                     bias_clo1=NA, bias_clo2=NA, rmse_clo1=NA, rmse_clo2=NA)
            
            
  
#x = 240        
            pb <- pbapply::timerProgressBar(max = nrow(out_x), style = 1, width = getOption("width")/4)
            progress <- function(x){setTimerProgressBar(pb, x)}
            opts <- list(progress = progress)
            outlist_x<-
              foreach(x = 1:nrow(out_x),
                      .packages = c("mice","plyr","tidyverse","data.table","dplyr","lpa.mi.src"), 
                      .inorder = TRUE, 
                      .options.snow = opts) %dopar% {
#  for(x in 1:nrow(out_x)){print(x)
                            
                         
                            # Get preliminaries
                            z_x = out_x$z[x]; pva_x = out_x$pva[x]; pm_x = out_x$pm[x]; type_x = out_x$data_type[x];
                            
                            
                            
                            # Get the parameters
                            parameters_x = parameters_combined_df %>% 
                              filter(rep==rep_x & z==z_x & pm==pm_x & data_type==type_x & pva==pva_x)
                            
                            return_x = out_x[x, ]
                            if(nrow(parameters_x)>0){

                                Qlist_x <- parameters_x %>% 
                                  select(paramHeader,param,LatentClass,est) %>% 
                                  Mplus2Qlist()
                                
                                Plist_x <- get_Plist(z_x, data_conditions)
                                
                                # Load the complete data
                                tmp1 = "complete"; tmp2=".RData"
                                load(paste0("H:/rdata-files/list-",tmp1," rep",rep_x," z",z_x,tmp2))
                                Ycomp_x = list_complete$dfcom %>% select(starts_with("Y"))
     
                                # Load the data
                                if(type_x=="Complete data"){tmp1 = "complete"; tmp2=".RData"}
                                if(type_x=="Observed data"){tmp1 = "observed"; tmp2=paste0(" pm",pm_x,".RData")}
                                if(type_x=="Imputation"){tmp1 = "imputed"; tmp2 = paste0(" pm",pm_x," pva",pva_x,".RData")}
                                load(paste0("H:/rdata-files/list-",tmp1," rep",rep_x," z",z_x,tmp2))
                                
                                P_cprob_x <- lpa.mi.src::cprobs(Y_i = Ycomp_x, pi_vec = as.numeric(Plist_x$pi[1,]), mu_mat = Plist_x$mu, S_array = Plist_x$S)
                                
                                if(type_x!="Imputation"){
                                  if(type_x=="Complete data"){Y_x = Ycomp_x}
                                  if(type_x=="Observed data"){Y_x = list_observed$list_obsdf$pm1 %>% select(starts_with("Y"))}
                                  Q_cprob_x <- lpa.mi.src::cprobs(Y_i = Y_x, pi_vec = Qlist_x$pi, mu_mat = Qlist_x$mu, S_array = Qlist_x$S)
                                  
                                  QmP = Q_cprob_x - P_cprob_x
                                  QmP_sq = QmP^2
                                  return_x$bias_cprob1 = mean(QmP[,1])
                                  return_x$bias_cprob2 = mean(QmP[,2])
                                  return_x$bias_cprob3 = mean(QmP[,3])
                                  return_x$rmse_cprob1 = mean(QmP_sq[,1])
                                  return_x$rmse_cprob2 = mean(QmP_sq[,2])
                                  return_x$rmse_cprob3 = mean(QmP_sq[,3])
                                  return_x$kl_cprob = sum(P_cprob_x*(log(Q_cprob_x) - log(P_cprob_x)))
                                  
                                  Q_lo_x = Q_cprob_x[,-3]/Q_cprob_x[,3]
                                  P_lo_x = P_cprob_x[,-3]/P_cprob_x[,3]
                                  lo_QmP = Q_lo_x - P_lo_x
                                  lo_QmP_sq = lo_QmP^2
                                  return_x$bias_clo1 = mean(lo_QmP[,1])
                                  return_x$bias_clo2 = mean(lo_QmP[,2])
                                  return_x$rmse_clo1 = mean(lo_QmP_sq[,1])
                                  return_x$rmse_clo2 = mean(lo_QmP_sq[,2])
                                  
                                  
                                } else {
                                  tmp_mids = list_imputed$obj_call[[pm_x]][[1]]
                                  list_Q_cprob_x<-
                                      lapply(X = 1:tmp_mids$m, 
                                             FUN = function(m){
                                               Y_m = mice::complete(tmp_mids, action = m) %>% select(starts_with("Y"));
                                               Q_cprob_m <- lpa.mi.src::cprobs(Y_i = Y_m, pi_vec = Qlist_x$pi, mu_mat = Qlist_x$mu, S_array = Qlist_x$S)
                                             })
                                  
                                  returnpp_x <- 
                                      lapply(X = 1:tmp_mids$m, 
                                             FUN = function(m){
                                               return_m = return_x; return_m$m = m; return_m$summary_type = "EAPP"
                                                Q_cprob_m <- list_Q_cprob_x[[m]]
                                                
                                                QmP = Q_cprob_m - P_cprob_x
                                                QmP_sq = QmP^2
                                                return_m$bias_cprob1 = mean(QmP[,1])
                                                return_m$bias_cprob2 = mean(QmP[,2])
                                                return_m$bias_cprob3 = mean(QmP[,3])
                                                return_m$rmse_cprob1 = mean(QmP_sq[,1])
                                                return_m$rmse_cprob2 = mean(QmP_sq[,2])
                                                return_m$rmse_cprob3 = mean(QmP_sq[,3])
                                                
                                                return_m$kl_cprob = sum(P_cprob_x*(log(Q_cprob_m) - log(P_cprob_x)))
                                                
                                                Q_lo_m = Q_cprob_m[,-3]/Q_cprob_m[,3]
                                                P_lo_m = P_cprob_x[,-3]/P_cprob_x[,3]
                                                lo_QmP = Q_lo_m - P_lo_m
                                                lo_QmP_sq = lo_QmP^2
                                                return_m$bias_clo1 = mean(lo_QmP[,1])
                                                return_m$bias_clo2 = mean(lo_QmP[,2])
                                                return_m$rmse_clo1 = mean(lo_QmP_sq[,1])
                                                return_m$rmse_clo2 = mean(lo_QmP_sq[,2])
                                                
                                                return(return_m)}
                                        ) %>% data.table::rbindlist() %>% data.frame()
                                    
                                  # Probability scale
                                    Q_cprob_x <- 0.*P_cprob_x
                                    for(m in 1:tmp_mids$m){Q_cprob_x<-list_Q_cprob_x[[m]]/tmp_mids$m + Q_cprob_x}
                                    QmP = Q_cprob_x - P_cprob_x
                                    QmP_sq = QmP^2
                                    return_x$bias_cprob1 = mean(QmP[,1])
                                    return_x$bias_cprob2 = mean(QmP[,2])
                                    return_x$bias_cprob3 = mean(QmP[,3])
                                    return_x$rmse_cprob1 = mean(QmP_sq[,1])
                                    return_x$rmse_cprob2 = mean(QmP_sq[,2])
                                    return_x$rmse_cprob3 = mean(QmP_sq[,3])
                                    return_x$kl_cprob = sum(P_cprob_x*(log(Q_cprob_x) - log(P_cprob_x)))
                                    
                                    Q_lo_x = Q_cprob_x[,-3]/Q_cprob_x[,3]
                                    P_lo_x = P_cprob_x[,-3]/P_cprob_x[,3]
                                    lo_QmP = Q_lo_x - P_lo_x
                                    lo_QmP_sq = lo_QmP^2
                                    return_x$bias_clo1 = mean(lo_QmP[,1])
                                    return_x$bias_clo2 = mean(lo_QmP[,2])
                                    return_x$rmse_clo1 = mean(lo_QmP_sq[,1])
                                    return_x$rmse_clo2 = mean(lo_QmP_sq[,2])
                                    return_x$summary_type = "pooled, prob. scale"
                                    
                                    
                                    # Odds ratio scale
                                    return_x2 = return_x
                                    Q_lo_x <- 0.*P_cprob_x[,-3]
                                    for(m in 1:tmp_mids$m){
                                      Q_lo_x<-log(list_Q_cprob_x[[m]][,-3]/cbind(list_Q_cprob_x[[m]][,3],list_Q_cprob_x[[m]][,3]))/tmp_mids$m  + Q_lo_x
                                    }
                                    Q_cprob_x = P_cprob_x + NA
                                    for(i in 1:nrow(Q_cprob_x)){
                                      Q_cprob_x[i,] = gamma2pi(Q_lo_x[i,])
                                    }
                                    QmP = Q_cprob_x - P_cprob_x
                                    QmP_sq = QmP^2
                                    return_x2$bias_cprob1 = mean(QmP[,1])
                                    return_x2$bias_cprob2 = mean(QmP[,2])
                                    return_x2$bias_cprob3 = mean(QmP[,3])
                                    return_x2$rmse_cprob1 = mean(QmP_sq[,1])
                                    return_x2$rmse_cprob2 = mean(QmP_sq[,2])
                                    return_x2$rmse_cprob3 = mean(QmP_sq[,3])
                                    return_x2$kl_cprob = sum(P_cprob_x*(log(Q_cprob_x) - log(P_cprob_x)))
                                    
                                    Q_lo_x = Q_cprob_x[,-3]/Q_cprob_x[,3]
                                    P_lo_x = P_cprob_x[,-3]/P_cprob_x[,3]
                                    lo_QmP = Q_lo_x - P_lo_x
                                    lo_QmP_sq = lo_QmP^2
                                    return_x2$bias_clo1 = mean(lo_QmP[,1])
                                    return_x2$bias_clo2 = mean(lo_QmP[,2])
                                    return_x2$rmse_clo1 = mean(lo_QmP_sq[,1])
                                    return_x2$rmse_clo2 = mean(lo_QmP_sq[,2])
                                    return_x2$summary_type = "pooled, log odds ratio scaled"
                                    
                                    return_x = rbind(return_x, return_x2)
                                    return_x = rbind(return_x, returnpp_x)
                                  
                                } #if(type_x)

                } #if(nrow(parameters_x)>0)
               return(return_x)
              }#end for(x = )
            
            out_x<-rbindlist(outlist_x) %>% data.frame()
            
            out_x$pva[out_x$data_type!="Imputation"] = NA
            save(out_x, file = paste0("H:/classify-accuracy-files/v3-classify-accuracy-rep",rep_x,".RData"))
            system(paste0('xcopy H:\\classify-accuracy-files\\v3-classify-accuracy-rep',rep_x,'.RData S:\\classify-accuracy-files /J /Y'))
  
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

