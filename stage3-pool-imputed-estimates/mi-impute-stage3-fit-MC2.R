rm(list = ls())

load("S:/environment-mi-impute-stage0 Jan 06 2019 13 29.RData")

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

system("mkdir E:\\pooled-files")


computer_name = "MC2"
pingpong_wd = paste0("S:/ping-pong")

Processors = 10
cl<-makeSOCKcluster(Processors)
doSNOW::registerDoSNOW(cl)    

for(rep in 1:Replications){
  
  if( !(paste0("rep",rep,".csv")%in%list.files(path = pingpong_wd)) ){
    
    tic = proc.time()
    write.csv(x = data.frame(computer = computer_name, total_time = NA), file = paste0(pingpong_wd,"/rep",rep,".csv"), row.names = FALSE)    
    
  
    system("mkdir H:\\est-files")
    system("mkdir H:\\tracker-files")
    system("mkdir H:\\pooled-files")
  
    system(paste0("xcopy S:\\tracker-files\\tracker-rep",rep,".RData H:\\tracker-files\\ /J /Y"))
    
    system(paste0("xcopy S:\\est-files\\rep",rep,".zip H:\\est-files\\ /J /Y"))
    system(paste0("mkdir H:\\est-files\\rep",rep))
    system(paste0("Bandizip.exe x -y -o:H:\\est-files\\rep",rep," H:\\est-files\\rep",rep,".zip"))
    
  # Load the proper tracker
    load(paste0("H:/tracker-files/tracker-rep",rep,".RData"))
   
    # Conduct the pooling
    tracker_df = tracker2_pool(tracker_df=tracker_df,
                pool_wd="H:/pooled-files",
                data_conditions=data_conditions,
                M_max = M,
                cl = cl,
                points.montecarlo = 1E5)
  
    # Save the resulting tracker and copy to S
    save(tracker_df, file = paste0("E:/pooled-tracker-rep",rep,".RData"))
    system(paste0("xcopy E:\\pooled-tracker-rep", rep, ".RData S:\\pooled-trackers /Y /J"))
    
    # Zip and the file to E and copy to S
    system(paste0("Bandizip.exe c -y E:\\pooled-files\\pooled-rep",rep,".zip H:\\pooled-files\\rep",rep))
    system(paste0("xcopy E:\\pooled-files\\pooled-rep",rep,".zip S:\\pooled-files /Y /J"))
    
    # Clean up 
    system("rm -r H:\\")
    rm(tracker_df)
    
    
    toc = proc.time()-tic; toc = round(toc[[3]],0)
    write.csv(x = data.frame(computer = computer_name, total_time = toc), file = paste0(pingpong_wd,"/rep",rep,".csv"), row.names = FALSE)    
  } #end if
}

stopCluster(cl)


