rm(list = ls())

library(tidyverse)
library(data.table)
library(pbapply)

results_wd = "S:/results"
cleaned_results_wd = "S:/cleaned-results"

files_vec = list.files(path = results_wd, pattern = ".RData")
param_files_vec = files_vec[startsWith(files_vec,"results-par")]
summary_files_vec = files_vec[startsWith(files_vec,"results-summaries")]

# Get summaries
clean_summary<-function(x){
  load(paste0(results_wd,"/",summary_files_vec[x]))
  cleaned_summary_df = summaries_df %>% select(Q,MM,ARIV,ARIU,Imputation_Procedure, 
                                             KL,se_KL, rep, z, pm, pva, kfit)
  return(cleaned_summary_df)
}
list_summaries = pblapply(X = 1:length(summary_files_vec), FUN=function(x){clean_summary(x)})
cleaned_summary_lpa_mi_df = rbindlist(list_summaries)


# Get parameters
clean_parameters<-function(x){
  load(paste0(results_wd,"/",param_files_vec[x]))
  cleaned_parameters_df = parameters_df %>% select(paramHeader,param,LatentClass,est,se,value,deviance,covered, 
                                                 rep,z,pm,pva,kfit,MM,data_type,Imputation_Procedure, data_type)
  return(cleaned_parameters_df)
}
list_parameters = pblapply(X = 1:length(param_files_vec), FUN = function(x){clean_parameters(x)})
cleaned_parameters_lpa_mi_df = rbindlist(list_parameters)

save(cleaned_summary_lpa_mi_df, file = "S:/cleaned-results/cleaned-summaries-lpa-mi-impute.RData")
save(cleaned_parameters_lpa_mi_df, file = "S:/cleaned-results/cleaned-parameters-lpa-mi-impute.RData")
