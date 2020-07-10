rm(list = ls())

library(plyr)
library(tidyverse)
library(pbapply)
library(data.table)
library(stringr)

load(file = "S:/environment-mi-impute-stage0 Jan 06 2019 13 29.RData")
data_conditions = transform(data_conditions, z = 1:nrow(data_conditions))
conditions_df = data.frame(z = 1:nrow(data_conditions)) %>%
  transform(Sample_Size = mapvalues(data_conditions$N,from = c(300,1200),c("Small N", "Large N"))) %>% 
  transform(Separation = mapvalues(data_conditions$MD, from = c(2.87,3.70), to = c("Weakly Sep.","Strongly Sep."))) %>% 
  transform(Mean_Differences = mapvalues(data_conditions$dX, from = c(T,F), to = c("Mean Diff.", "No Mean Diff."))) %>% 
  transform(Association = mapvalues(data_conditions$rho_YX, from = c(0.4,0), to = c("Correlation", "No Correlation"))) %>% 
  transform(Mixing = mapvalues(data_conditions$class_size, from = c("A","B"), to = c("Equal Mixing","Unequal Mixing"))) %>% 
  transform(Moderation = mapvalues(data_conditions$C_modifies_YX, from = c(T,F),to = c("Moderation","No Moderation")))
            


trackers_wd = "H:/pooled-trackers"
poolfiles_wd = "H:/pooled-files"
results_wd = "H:/results"


tmp<-pblapply(X = 1:Replications, 
              FUN = function(x){
                  load(paste0(trackers_wd,"/pooled-tracker-rep",x,".RData"))
                  # delete the estwd.y and estfolder.y
                  cols_del = which(names(tracker_df)=="estwd.y" | names(tracker_df)=="estfolder.y")
                  tracker_df = tracker_df[,-cols_del]
                  names(tracker_df)[which(names(tracker_df)=="estwd.x")] = "estwd"
                  names(tracker_df)[which(names(tracker_df)=="estfolder.x")] = "estfolder"
                  
                  # Generate a pick1 variable
                  pick1_df = ddply(subset(tracker_df, data_type == "Imputed data" & converged == TRUE & M_sufficient==TRUE), .(rep,z,pm,pva,kfit), summarize,
                                   m = min(m)) %>% 
                             transform(pick1 = 1) 
                  tracker_df = left_join(x = tracker_df, y = pick1_df)
                  
                  # Update the name of the poolfolder 
                  tracker_df = transform(tracker_df, poolfolder = str_replace(string = poolfolder, pattern = "/rep","pooled-rep"))
                  
                  
                  # rows of pick1 
                  tmp_df = tracker_df %>% filter(pick1==1) %>% select(rep,z,pm,pva,kfit,pick1,poolwd,poolfolder,poolfile) %>% 
                     transform(Q=NA,MM=NA,ARIV=NA,ARIU=NA,KL=NA,se_KL=NA)
                  for(i in seq(1,nrow(tmp_df))){
                    load(with(tmp_df[i,], paste(poolwd,poolfolder,poolfile,sep = "/")))
                    
                    tmp_df$Q[i] = pooled_list$additional$Q
                    tmp_df$MM[i] = pooled_list$additional$MM
                    tmp_df$ARIV[i] = pooled_list$additional$ARIV
                    tmp_df$ARIU[i] = pooled_list$additional$ARIU
                    tmp_df$KL[i] = pooled_list$additional$KL
                    tmp_df$se_KL[i] = pooled_list$additional$se_KL
                    
                    rm(pooled_list)
                  }
                  
                  summaries_df = left_join(x = tmp_df, y = conditions_df)
                  summaries_df = transform(summaries_df, 
                                           Imputation_Procedure = mapvalues(summaries_df$pva, from = 1:4, to = c("Amelia","PMM","CART","RF")))
                  save(summaries_df, file = paste0(results_wd,"/results-summaries-rep",x,".RData"))
                  
                  tmp_list = lapply(X = seq(1,nrow(summaries_df)), 
                         FUN =function(i){
                                load(with(summaries_df[i,], paste(poolwd,poolfolder,poolfile,sep = "/")))
                                parameters_df = pooled_list$parameters %>%
                                                transform(rep = summaries_df$rep[i], 
                                                          z =  summaries_df$z[i], 
                                                          pm =  summaries_df$pm[i], 
                                                          pva =  summaries_df$pva[i], 
                                                          kfit =  summaries_df$kfit[i], 
                                                          MM = summaries_df$MM[i],
                                                          Sample_Size = summaries_df$Sample_Size[i], 
                                                          Separation = summaries_df$Separation[i], 
                                                          Mean_Differences = summaries_df$Mean_Differences[i], 
                                                          Association = summaries_df$Association[i], 
                                                          Mixing = summaries_df$Mixing[i], 
                                                          Moderation = summaries_df$Moderation[i],
                                                          Imputation_Procedure = summaries_df$Imputation_Procedure[i])
                    }) 
                  parameters_df = rbindlist(tmp_list)
                  save(parameters_df, file = paste0(results_wd,"/results-parameters-rep",x,".RData"))
        
    return(NULL)
})

