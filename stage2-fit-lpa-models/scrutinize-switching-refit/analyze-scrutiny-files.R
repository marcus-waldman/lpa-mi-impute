rm(list = ls())

require(pbapply)
require(plyr)
require(tidyverse)
require(data.table)

# Load in the scrutiny file
files_scrutiny = list.files("S:/scrutiny-files/", full.names = TRUE)


list_scrutiny<-pblapply(X = files_scrutiny, 
                        FUN = function(x){load(x); return(scrutinize_df)})
scrutiny_df = rbindlist(list_scrutiny)
rm(list_scrutiny)
gc()


# Load in the tracker files
files_tracker = list.files("S:/tracker-files/", full.names = TRUE)
list_trackers = pblapply(X = files_tracker, 
                        FUN = function(x){load(x);return(tracker_df)})
tracker_df = rbindlist(list_trackers)
rm(list_trackers)
gc()


# Summarize the scrutiny file
head(scrutiny_df)
scrutiny_df %>% transform(max_L1_ifom_est = max(c(L1_mu_ifom_est, L1_var_ifom_est))) %>%
                transform(max_L1_inp_ifom = max(c(L1_mu_inp_ifom, L1_var_inp_ifom))) %>% 
                transform(ratio_Rcond = Rcond_est/Rcond_tracker)


