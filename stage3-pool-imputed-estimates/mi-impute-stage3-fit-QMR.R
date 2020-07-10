rm(list = ls())

load("A:/environment-mi-impute-stage0 Jan 06 2019 13 29.RData")

library(tidyverse)
library(lpa.mi.src)
library(data.table)
library(snow)
library(doSNOW)
library(foreach)


# Construct the in the methods-list
M = 100
VerboseImpute = FALSE
methods_list = list(procedure = list(NULL), name = list(NULL), args = list(NULL))
pva_vec = c(1:4)
methods_list$procedure[[1]] = "amelia";         methods_list$name[[1]] = "Amelia";       methods_list$args[[1]] = list(m = M, p2s = as.numeric(VerboseImpute), tolerance = 1E-4)
methods_list$procedure[[2]] = "mice";           methods_list$name[[2]] = "PMM";          methods_list$args[[2]] = list(method = "pmm", maxit = 50, m = M, printFlag = VerboseImpute)
methods_list$procedure[[3]] = "mice";           methods_list$name[[3]] = "CART";         methods_list$args[[3]] = list(method = "cart", maxit = 50, m = M, printFlag = VerboseImpute)
methods_list$procedure[[4]] = "mice";           methods_list$name[[4]] = "RF";           methods_list$args[[4]] = list(method = "rf", maxit = 25, ntree = 10, m = M, printFlag = VerboseImpute)

est_zip_wd = "A:/est-files"
pool_wd = "A:/pool-files"


Processors = 14
cl<-makeSOCKcluster(Processors)
doSNOW::registerDoSNOW(cl)    


for(rep in 5:Replications){

# Load the proper tracker
  load(paste0("A:/tracker-files/tracker-rep",rep,".RData"))
  
  # print(paste0("Unzipping for replication ", rep))
  # # Unzip the file
   est_unzip_wd = "H:/est-files/"
  # unzip(zipfile = paste0(est_zip_wd,"/rep",rep,".zip"), 
  #       exdir = paste0("H:/est-files/rep",rep))
  # 
  # Conduct the pooling
  tracker_df = tracker2_pool(tracker_df=tracker_df,
              pool_wd=pool_wd,
              est_zip_wd=est_zip_wd,
              est_unzip_wd=est_unzip_wd,
              data_conditions=data_conditions,
              M_max = M,
              cl = cl,
              points.montecarlo = 1E5)

  # Save the resulting tracker
  save(tracker_df, file = paste0("H:/tracker-files/tracker-rep",rep,".RData"))

# Clean up
rm(tracker_df)
  
}

stopCluster(cl)


