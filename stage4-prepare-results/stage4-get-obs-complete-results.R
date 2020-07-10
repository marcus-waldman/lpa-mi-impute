rm(list = ls())

library("stringr")
library("plyr")
library("tidyverse")
library("pbapply")
library("data.table")
library("lpa.mi.src")

load("F:/QMR-Recovery/environment-mi-impute-stage0 Jan 06 2019 13 29.RData")

wd = getwd()

for (x in 1:500){

    # Est-files copy over and extract
    system('mkdir E:\\est-files')
    system(paste0('mkdir E:\\est-files\\rep', x))
    system(paste0('xcopy F:\\SRS_24-02-2019_at_00.32.23\\diss-drive\\est-files\\rep', x, '.zip E:\\est-files\\rep',x,' /J /Y'))
    system(paste0('Bandizip.exe x -y E:\\est-files\\rep', x,'\\rep',x,'.zip'))
    system(paste0('rm E:\\est-files\\rep',x,'\\rep',x,'.zip'))
    
    
    # Tracker file copy over and extract
    system('mkdir E:\\tracker-files')
    system(paste0('xcopy F:\\QMR-Recovery\\pooled-trackers\\pooled-tracker-rep', x,'.RData E:\\tracker-files /J /Y'))
    load(paste0("E:\\tracker-files\\pooled-tracker-rep",x,".RData"))
    tracker_df = subset(tracker_df, data_type != "Imputed data" & converged == TRUE)
    tracker_df = transform(tracker_df, estwd.x = str_replace(string = estwd.x, pattern = "H:/", replacement = "E:/"))
    
    
    # Results for parameters
    system('mkdir E:\\results')
    system(paste0('xcopy F:\\QMR-Recovery\\results\\results-parameters-rep', x,'.RData E:\\results /J /Y'))
    load(paste0("E:\\results\\results-parameters-rep",x,".RData"))
    parameters_df = transform(parameters_df, data_type = "Imputation")
    
    pick1_df = ddply(.data = parameters_df, .(rep,z,pm,kfit,Sample_Size,Separation,Mean_Differences,Association,Mixing,Moderation,LatentClass,paramHeader,param), summarize, 
                     value = mean(value))
    
    tmp_df = pblapply(X = seq(1,nrow(tracker_df)), 
                      FUN =function(i){
                        load(with(tracker_df[i,], paste(estwd.x,estfolder.x,estfile,sep = "/")))
                        tmp2_df = list_estimates$out_Mplus$parameters$unstandardized %>% 
                                 transform(rep = tracker_df$rep[i], 
                                           z = tracker_df$z[i], 
                                           kfit = tracker_df$kfit[i],
                                           data_type = tracker_df$data_type[i]) 
                        tmp2_df = subset(tmp2_df, !endsWith(tmp2_df$paramHeader, ".WITH"))
                        
                        return(tmp2_df)
                      }) %>% 
                rbindlist() %>% 
                left_join(y = pick1_df, by = c("paramHeader", "param", "LatentClass", "rep", "z", "kfit")) %>%
                transform(deviance = est - value,
                          covered = ifelse(value>=(est-1.959964*se) & value<(est+1.959964*se),TRUE,FALSE), 
                          pva = NA, 
                          MM = NA, 
                          Imputation_Procedure = ifelse(data_type == "Complete data", "No Missingness", "FIML"))
    tmp_df = tmp_df[,names(parameters_df)]
    parameters_df = rbind(tmp_df, parameters_df)
    save(parameters_df, file = paste0("E:/results/results-parameters-rep",x,".RData"))
    rm(tmp_df)

    
    # Results for KL divergence
    system(paste0('xcopy F:\\QMR-Recovery\\results\\results-summaries-rep', x,'.RData E:\\results /J /Y'))
    load(paste0("E:\\results\\results-summaries-rep",x,".RData"))
    
    pick1_df = ddply(.data = parameters_df, .(z,pm,kfit,Sample_Size,Separation,Mean_Differences,Association,Mixing,Moderation), summarize, 
                     rep = mean(rep))
    
    
    tmp_df = pblapply(X = seq(1,nrow(tracker_df)), 
                        FUN = function(i){
                          load(with(tracker_df[i,], paste(estwd.x,estfolder.x,estfile,sep = "/")))
                          Plist_x = lpa.mi.src::get_Plist(z = tracker_df$z[i], data = data_conditions)
                          Qlist_x = lpa.mi.src::Mplus2Qlist(params_df = list_estimates$out_Mplus$parameters$unstandardized)
                          KL_df = lpa.mi.src::KLmixmvrnorm(P = Plist_x, Q = Qlist_x)
                          tmp2_df = data.frame(rep = tracker_df$rep[i], 
                                              z = tracker_df$z[i], 
                                              kfit = tracker_df$kfit[i],
                                              data_type = tracker_df$data_type[i],
                                              KL = KL_df$est, 
                                              se_KL = KL_df$se)
                          
                          return(tmp2_df)
                        }) %>% 
                rbindlist() %>% 
                left_join(y = pick1_df, by = c("rep", "z", "kfit")) %>% 
                transform(ARIV = NA, ARIU = NA, pva = NA, pick1 = NA, poolwd = NA, poolfolder = NA, poolfile = NA, Q = NA, MM = NA, 
                          Imputation_Procedure = ifelse(data_type == "Complete data", "No Missingness", "FIML"))
    tmp_df = tmp_df[,names(summaries_df)]
    summaries_df = rbind(tmp_df, summaries_df)
    rm(tmp_df)
    save(summaries_df, file = paste0("E:/results/results-summaries-rep",x,".RData"))
    
    
    # Copy over the est-files
    system(paste0('Bandizip.exe c -y E:\\est-files\\est-files-rep',x,'.zip E:\\est-files\\rep',x,'\\'))
    system(paste0('xcopy E:\\est-files\\est-files-rep',x,'.zip F:\\diss-drive\\est-files /Y /J'))
    
    # Copy over the results
    system(paste0('xcopy E:\\results\\results-parameters-rep',x,'.RData F:\\diss-drive\\results\\ /Y /J'))
    system(paste0('xcopy E:\\results\\results-summaries-rep',x,'.RData F:\\diss-drive\\results\\ /Y /J'))
    
    # Clean up
    rm(pick1_df)
    rm(tracker_df)
    rm(parameters_df)
    rm(summaries_df)
    
    system(paste0('rm -r E:\\'))
    
}



