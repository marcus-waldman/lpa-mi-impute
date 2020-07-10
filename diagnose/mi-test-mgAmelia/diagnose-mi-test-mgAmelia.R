rm(list = ls())

library(plyr)
library(tidyverse)
library(lpa.mi.src)

setwd("C:/Users/marcu/Dropbox/Dissertation/saved/data-gen/test-mgAmelia/test-mgAmelia Jan 03 2019 10 11")
load("environment.RData")


z = 1
rep = 1
pop_params_z = get_FMM_params(z = z, data_conditions = data_conditions)


# Unzip and load complete and observed data
file_unzip = paste0("list-stage0 rep", rep, " z",z,".RData")
dir_unzip = paste0(temp_wd,"/stage0")
unzip(zipfile = paste(stage0_dir, stage0_zip, sep = "/"), 
      files = file_unzip, 
      exdir = dir_unzip)
load(paste0(dir_unzip,"/",file_unzip))



list_get_obs = list_observed
list_get_complete = list_complete
save_it = TRUE
rep = rep 
temp_wd_rep_vec = datfiles_wd_rep_vec


# Generate imputations
list_imputed = lpa.mi.src::get_imputed_data(z = z,
                                            list_get_obs = list_observed,
                                            list_get_complete = list_complete,
                                            methods_list = methods_list,
                                            data_conditions = data_conditions,
                                            pctmiss_vec = pctmiss_vec,
                                            save_it = FALSE)
#              print("  >> list_imputed constructed")
