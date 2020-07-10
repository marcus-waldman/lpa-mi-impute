rm(list = ls())

library(plyr)
library(tidyverse)
library(data.table)
library(pbapply)

pattmix_results_wd = "S:/lpa-pattmix-files/results"
impute_results_wd = "S:/cleaned-results"
dropbox_stage4c_wd = "D:/Dropbox/Dissertation/lpa-mi-impute/stage4c-combine-results"

setwd(impute_results_wd)
load(file = "cleaned-parameters-lpa-mi-impute.RData") 
load(file = "cleaned-summaries-lpa-mi-impute.RData")

setwd(pattmix_results_wd)
load(file = "results-parameters-lpa-pattmix.RData")
load(file = "results-summary-lpa-pattmix.RData")


# Apppend the summaries
x = names(cleaned_summary_lpa_mi_df); y = names(summaries_lpa_pattmix_df)
dset_summary = union(setdiff(x,y), setdiff(y,x))
dset_summary[!(dset_summary %in% names(cleaned_summary_lpa_mi_df))] #Variables to add to cleaned_summary.
cleaned_summary_lpa_mi_df = cleaned_summary_lpa_mi_df %>% transform(points_montecarlo = NA, data_type = "NA")

dset_summary[!(dset_summary %in% names(summaries_lpa_pattmix_df))] #Variables to add to summaries_lpa
summaries_lpa_pattmix_df = transform(summaries_lpa_pattmix_df, 
                                     Q = NA, MM = NA, ARIV = NA, ARIU = NA, Imputation_Procedure = "Pattern Mixture", 
                                     pva = NA)
setdiff(names(summaries_lpa_pattmix_df), names(cleaned_summary_lpa_mi_df))
setdiff(names(cleaned_summary_lpa_mi_df), names(summaries_lpa_pattmix_df))

class(cleaned_summary_lpa_mi_df)
summaries_lpa_pattmix_df = summaries_lpa_pattmix_df%>% as.data.table()
class(summaries_lpa_pattmix_df)

summaries_combined_df = bind_rows(cleaned_summary_lpa_mi_df,summaries_lpa_pattmix_df) 
rm(summaries_lpa_pattmix_df); rm(cleaned_summary_lpa_mi_df);
save(summaries_combined_df, file = paste0(dropbox_stage4c_wd,"/summaries-combined-results-lpa-mi-impute.RData"))


# Append the Parameters
x = names(parameters_lpa_pattmix_df); y = names(cleaned_parameters_lpa_mi_df)
dset_summary = union(setdiff(x,y), setdiff(y,x))

dset_summary[!(dset_summary %in% names(cleaned_parameters_lpa_mi_df))] #Variables to add to cleaned_parameters.

dset_summary[!(dset_summary %in% names(parameters_lpa_pattmix_df))] #Variables to add to cleaned_parameters.
parameters_lpa_pattmix_df = parameters_lpa_pattmix_df %>% 
  transform(se = NA, covered = NA, pva = NA, MM = NA, Imputation_Procedure = "Pattern Mixture")

class(parameters_lpa_pattmix_df)
parameters_lpa_pattmix_df = as.data.table(parameters_lpa_pattmix_df)
class(cleaned_parameters_lpa_mi_df)

parameters_combined_df = bind_rows(cleaned_parameters_lpa_mi_df, 
                                   parameters_lpa_pattmix_df)
parameters_combined_df = parameters_combined_df %>% filter(!endsWith(paramHeader,".WITH"))
save(parameters_combined_df, file = paste0(dropbox_stage4c_wd,"/parameters-combined-results-lpa-mi-impute.RData"))

