#  .rs.restartR()

rm(list = ls())

library(plyr)
library(tidyverse)
library(lpa.mi.src)
library(snow)
library(doSNOW)
library(foreach)
library(data.table)
library(pbapply)

computer_name = "MC3"
pingpong_wd = paste0("S:/ping-pong")

Processors = 10
cl<-makeSOCKcluster(Processors)
doSNOW::registerDoSNOW(cl) 



# Number of classes fit 
K_vec = 3;

# Set the working directory
setwd("H:/")

# Load the working environment
load("S:/environment-mi-impute-stage0 Jan 06 2019 13 29.RData")
#rm(methods_list)


z_vec = 1:40
pm_vec = 1
pva_vec = 1:4
M = 100
dk_vec = 0
starts0 = 20

repseq = seq(250,500,by=1)

#rep = repseq[2]
for (rep in repseq){
  if( !(paste0("rep",rep,".csv")%in%list.files(path = pingpong_wd)) ){
    
    tic = proc.time()
    write.csv(x = data.frame(computer = computer_name, total_time = NA), file = paste0(pingpong_wd,"/rep",rep,".csv"), row.names = FALSE)    
    
    system("mkdir H:\\est-files")
    system("mkdir H:\\out-files")
    system("mkdir H:\\dat-files")
    
    # Copy needed files from S drive
    system(paste0('xcopy S:\\dat-files\\dat-files-rep',rep,'.zip E:\\ /Y /J'))
    system(paste0('xcopy S:\\tracker-files\\tracker-rep',rep,'.RData E:\\ /Y'))
    
    # Extract dat files H drive (working directory)
    system(paste0('Bandizip x -y -o:H:\\dat-files E:\\dat-files-rep',rep,'.zip'))
    
    # Load the tracker and the subset to be refit, along with only the relavant columns 
    load(paste0("E:/tracker-rep",rep,".RData"))
    tracker_df = subset(tracker_df, pass_scrutiny == FALSE)[,1:20]
    
    # Clean the tracker
    library(stringr)
    tracker_df = transform(tracker_df, datwd = str_replace(string = datwd, pattern = "D:/",replacement = "H:/"))
    tracker_df = transform(tracker_df, outwd = str_replace(string = outwd, pattern = "D:/",replacement = "H:/"))
    
    
    # Create out-files folders
    print("Creating nessary out-files directories:")
    tmp<-pblapply(X = 1:nrow(tracker_df), 
                  FUN = function(x){
                    track_x = tracker_df[x,];
                    with(track_x,
                    dir.create(path = paste0(outwd,outfolder), 
                               recursive = T, 
                               showWarnings = F)
                         );
                  }
    )
    
    tracker_df = tracker2_fit_til_convergence(data_conditions = data_conditions,
                                              tracker_df = tracker_df,
                                              cl = cl)

    tracker_df$converged = FALSE
    tracker_df$converged[tracker_df$normal==TRUE & tracker_df$poor_Rcond==FALSE] = TRUE
    
    tracker_df = tracker2_switch(tracker_df = tracker_df,
                                 data_conditions = data_conditions, 
                                 cl = cl)
    
    tracker_df = tracker2_save_readMplus(tracker_df=tracker_df, estwd = "H:/est-files", cl = cl)
    
    # Update the est-files
    system(paste0('xcopy S:\\est-files\\rep',rep,'.zip E:\\est-files\\ /Y /J'))
    system(paste0('Bandizip.exe a -r -aoa -y E:\\est-files\\rep', rep,'.zip H:\\est-files\\rep',rep,'\\*.RData'))
    system(paste0('xcopy E:\\est-files\\rep', rep,'.zip S:\\est-files /Y /J'))
    system(paste0('rm E:\\est-files\\rep',rep,'.zip'))
    system('rm -r H:\\est-files')
    
    # Update the out-files
    system(paste0('xcopy S:\\out-files\\out-files-rep', rep, ".zip E:\\out-files\\ /Y /J"))
    system(paste0('Bandizip.exe a -r -aoa -y -l:9 E:\\out-files\\out-files-rep', rep,'.zip H:\\out-files\\rep',rep,'\\*.inp'))
    system(paste0('Bandizip.exe a -r -aoa -y -l:9 E:\\out-files\\out-files-rep', rep,'.zip H:\\out-files\\rep',rep,'\\*.out'))
    system(paste0('xcopy E:\\out-files\\out-files-rep', rep,'.zip S:\\out-files /Y /J'))
    system(paste0('rm E:\\out-files\\out-files-rep',rep,'.zip'))
    system('rm -r H:\\out-files')

    
    # Send over the tracker-file
    save(tracker_df, file = paste0("E:/refits-tracker-rep",rep,".RData"))
    save(tracker_df, file = paste0("S:/tracker-files/refits/refits-tracker-rep",rep,".RData"))
    
    # Clean up the dat-files
    system(paste0('rm E:\\dat-files-rep',rep,'.zip'))
    system('rm -r H:\\dat-files')

    
    
    toc = proc.time()-tic; toc = round(toc[[3]],0)
    write.csv(x = data.frame(computer = computer_name, total_time = toc), file = paste0(pingpong_wd,"/rep",rep,".csv"), row.names = FALSE)    
    
  } #end if
} #end for rep in 


stopCluster(cl)