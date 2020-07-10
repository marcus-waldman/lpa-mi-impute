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

computer_name = "MC3"
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
system("rm -r H:\\classify-accuracy-files")

Processors = 10
cl<-makeSOCKcluster(Processors)
doSNOW::registerDoSNOW(cl)   


# Load in the results
setwd(results_wd)
load(file ="parameters-combined-results-lpa-mi-impute.RData")
parameters_combined_df$pva[parameters_combined_df$data_type=="Complete data"] = -1
parameters_combined_df$pva[parameters_combined_df$data_type=="Observed data"] = 0

for(rep_x in sample(1:500,500,replace=F)){

  if( !(paste0("rep",rep_x,".csv")%in%list.files(path = pingpong_wd)) ){
    
    print(paste0("Replication: ", rep_x))
    tic = proc.time()
    write.csv(x = data.frame(computer = computer_name, total_time = NA), file = paste0(pingpong_wd,"/rep",rep_x,".csv"), row.names = FALSE)    
    
  
          # Make replication directories
          system("rm -r H:\\rdata-files")
          system('mkdir H:\\rdata-files')
          
          system("rm-r H:\\classify-accuracy-files")
          system('mkdir H:\\classify-accuracy-files')
          
          # Copy over the complete data files
          system(paste0('xcopy "S:\\rdata-files\\list-complete rep',rep_x,' *.RData" H:\\rdata-files'))
          
          # copy over the observed data files
          system(paste0('xcopy "S:\\rdata-files\\list-observed rep',rep_x,' *.RData" H:\\rdata-files'))
          
          # copy over the imputed data files
          system(paste0('xcopy "S:\\rdata-files\\list-imputed rep',rep_x,' *.RData" H:\\rdata-files'))
          
          
          # create a subpopulation list
          list_subpop<-
            lapply(X = z_vec, FUN = function(zz){
              tmp1 = "complete"; tmp2=".RData";
              load(paste0("H:/rdata-files/list-",tmp1," rep",rep_x," z",zz,tmp2))
              return(list_complete$dfcom %>% select("subpop"))
            })
          
          
          
          out_x = expand.grid(z = z_vec, pva = c(-1:4), pm=1) %>% data.frame() %>% transform(data_type = NA)
          out_x$data_type[out_x$pva==-1] = "Complete data"
          out_x$data_type[out_x$pva==0] = "Observed data"
          out_x$data_type[out_x$pva>0] = "Imputation"
          out_x = out_x %>% transform(kappa1c1=NA, kappa1c2=NA, kappa1c3=NA, kappa2c1=NA, kappa2c2=NA, kappa2c3=NA, kappa3c1=NA, kappa3c2=NA, kappa3c3=NA)
          
          
          pb <- pbapply::timerProgressBar(max = nrow(out_x), style = 1, width = getOption("width")/4)
          progress <- function(x){setTimerProgressBar(pb, x)}
          opts <- list(progress = progress)
          outlist_x<-
            foreach(x = 1:nrow(out_x),
                    .packages = c("mice","plyr","tidyverse","data.table","dplyr","lpa.mi.src"), 
                    .inorder = TRUE, 
                    .options.snow = opts) %dopar% {
          #for(x in 1:nrow(out_x)){print(x)
                          z_x = out_x$z[x]; pva_x = out_x$pva[x]; pm_x = out_x$pm[x]; type_x = out_x$data_type[x]; 
                          
                          # Get the parameters
                          parameters_x = parameters_combined_df %>% 
                            filter(rep==rep_x & z==z_x & pm==pm_x & data_type==type_x & pva==pva_x)
                          
                          if(nrow(parameters_x)>0){
                          
                          
                              Qlist_x <- parameters_x %>% 
                                select(paramHeader,param,LatentClass,est) %>% 
                                Mplus2Qlist()
                              
                              
                              # Load the data
                              if(type_x=="Complete data"){tmp1 = "complete"; tmp2=".RData"}
                              if(type_x=="Observed data"){tmp1 = "observed"; tmp2=paste0(" pm",pm_x,".RData")}
                              if(type_x=="Imputation"){tmp1 = "imputed"; tmp2 = paste0(" pm",pm_x," pva",pva_x,".RData")}
                              load(paste0("H:/rdata-files/list-",tmp1," rep",rep_x," z",z_x,tmp2))
                              
                              if(type_x!="Imputation"){
                                if(type_x=="Complete data"){Y_x = list_complete$dfcom %>% select(starts_with("Y"))}
                                if(type_x=="Observed data"){Y_x = list_observed$list_obsdf$pm1 %>% select(starts_with("Y"))}
                                cprob_x <- lpa.mi.src::cprobs(Y_i = Y_x, pi_vec = Qlist_x$pi, mu_mat = Qlist_x$mu, S_array = Qlist_x$S)
                              
                              } else {
                                tmp_mids = list_imputed$obj_call[[pm_x]][[1]]
                                tmp_cprobs<-
                                  lapply(X = 1:tmp_mids$m, 
                                         FUN = function(m){
                                           Y_x = mice::complete(tmp_mids, action = m) %>% select(starts_with("Y"));
                                           cprob_x = lpa.mi.src::cprobs(Y_i = Y_x, pi_vec = Qlist_x$pi, mu_mat = Qlist_x$mu, S_array = Qlist_x$S) %>% data.frame() %>% transform(id = 1:nrow(Y_x), m = m)
                                           return(cprob_x)
                                         } ) %>% data.table::rbindlist()
                                cprob_x = tmp_cprobs %>% group_by(id) %>% summarise(X1 = mean(X1), X2 = mean(X2), X3 = mean(X3)) %>% select(X1,X2,X3) %>% data.frame()
                              }
                              modal_x = apply(cprob_x, 1, which.max)
                              table_x = table(list_subpop[[z_x]]$subpop, modal_x)
                              out_x$kappa1c1[x] = table_x[1,1]
                              out_x$kappa1c2[x] = table_x[1,2]
                              out_x$kappa1c3[x] = table_x[1,3]
                              out_x$kappa2c1[x] = table_x[2,1]
                              out_x$kappa2c2[x] = table_x[2,2]
                              out_x$kappa2c3[x] = table_x[2,3]
                              out_x$kappa3c1[x] = table_x[3,1]
                              out_x$kappa3c2[x] = table_x[3,2]
                              out_x$kappa3c3[x] = table_x[3,3]
                              
                              return(out_x[x, ])
              
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
          system("rm-r H:\\classify-accuracy-files")

  }#end if(ping)
}# end for rep=

stopCluster(cl)

