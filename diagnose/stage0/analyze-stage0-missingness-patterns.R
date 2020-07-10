rm(list = ls())

library(plyr)
library(tidyverse)
library(lpa.mi.src)

stage0_wd = "C:/Users/marcu/Dropbox/Dissertation/data-manipulation/lists-stage0"
zip_stage0 = "Rdata-mi-impute-stage0 Jan 05 2019 17 31.zip"

tmpfolder = ("C:/Users/marcu/Documents/zip-files/super-tmp")

setwd(stage0_wd)
unzip(zipfile = zip_stage0, exdir = tmpfolder)

setwd(tmpfolder)
files_vec=list.files(path = tmpfolder, pattern = ".RData")
load(files_vec[startsWith(files_vec,"environment")])
files_vec = files_vec[startsWith(files_vec, "list-stage0")]

repz_df = data.frame(expand.grid(z = 1:40, rep = 1:500))

getot<-function(x){
  rep = repz_df$rep[x]
  z = repz_df$z[x]
  load(files_vec[endsWith(files_vec, paste0("rep",rep," z",z,".RData"))])
  #obs_tmp = data.frame(obs_rate = min(list_observed$obs_rates[,-1]), obs_counts = min(list_observed$obs_counts[,-1]), 
  #                     rep = rep, z = z)
  #obs_tmp = list_observed$obs_rates %>% transform(z = z, rep = rep)
  obs_tmp = list_observed$obs_counts %>% transform(z = z, rep = rep)
  return(obs_tmp)
}
tmp_list = lapply(X = 1:nrow(repz_df), FUN = "getot")

library(data.table)

df = rbindlist(tmp_list)

data_conditions = transform(data_conditions, z = 1:Ztot)
nrow(df)
df = merge(df, data_conditions)
nrow(df)



ddply(df, .(z), summarize,
      Y1 = min(Y1), Y2 = min(Y2), Y3 = min(Y3), Y4 = min(Y4))

g