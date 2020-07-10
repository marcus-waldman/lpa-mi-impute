rm(list = ls())

#
computer_name = "MC1"
Processors = 10


# Copy over the environment variable
system('xcopy "S:\\environment-mi-impute-stage0 Jan 06 2019 13 29.RData" H:\\ /Y')


# Make the necessary directories
system('mkdir H:\\est-files')
system('mkdir H:\\out-files')
system('mkdir H:\\tracker-files')
system('mkdir H:\\scrutiny-files')

load("H:/environment-mi-impute-stage0 Jan 06 2019 13 29.RData")

library(stringr)
require(plyr)
require(tidyverse)
require(pbapply)
require(lpa.mi.src)
require(data.table)
require(snow)
require(doSNOW)
require(foreach)

syntax2Plist<-function(syntax_x){
  library(stringr)
  library(lpa.mi.src)
  
  syntax_x = toupper(syntax_x)
  th = c(which(str_detect(string = syntax_x, pattern = c("%C#"))),
         length(syntax_x)+1)
  K_x = length(th)-1
  template_list = list(class = NA, syntax = NA, cov = NA, mean = NA, var = NA)
  syntax_list<-lapply(X = 1:K_x, FUN = function(k){
    list_k = template_list;
    syntax_k = toupper(str_trim(syntax_x[seq(th[k], th[k+1]-1)]));
    
    idx_with_k = which(str_detect(string = syntax_k, pattern = "WITH"));
    idx_mean_k = which(str_detect(string = syntax_k, pattern = "\\["));
    idx_var_k = setdiff(which(str_detect(string = syntax_k, pattern = "Y")), 
                        c(idx_with_k, idx_mean_k));
    
    list_k$class = k
    list_k$syntax = syntax_k;
    list_k$cov = syntax_k[idx_with_k]; 
    list_k$mean = syntax_k[idx_mean_k];
    list_k$var = syntax_k[idx_var_k]
    
    return(list_k);
  })
  
  
  
  Plist_x = get_Plist(1, data_conditions = data_conditions)
  Plist_x$mu = NA*Plist_x$mu 
  Plist_x$S = NA*Plist_x$S 
  Plist_x$pi = NA*Plist_x$pi 
  
  for (k in 1:K_x){
    
    list_k = syntax_list[[k]]
    
    # means
    tmp = list_k$mean
    tmp = str_remove(string = tmp, pattern = c("\\["))
    tmp = str_remove(string = tmp, pattern = c("\\];"))
    for (j in 1:length(tmp)){
      tmp = str_remove(tmp, paste0("Y",j))
    }
    tmp = str_remove(string = tmp, pattern = c("\\*"))
    tmp = str_remove(string = tmp, pattern = c("\\@"))
    tmp = str_trim(tmp)
    Plist_x$mu[,k] = as.numeric(tmp)
    
    
    # cov
    tmp = list_k$cov
    tmp = str_remove(string = tmp, pattern = c("WITH"))
    tmp = str_remove(string = tmp, pattern = c(";"))
    for (j in 1:length(tmp)){
      tmp = str_remove(tmp, paste0("Y",j))
    }
    tmp = str_remove(string = tmp, pattern = c("\\*"))
    tmp = str_remove(string = tmp, pattern = c("\\@"))
    tmp = str_trim(tmp)
    
    # var
    tmp = list_k$var
    tmp = str_remove(string = tmp, pattern = c(";"))
    for (j in 1:length(tmp)){
      tmp = str_remove(tmp, paste0("Y",j))
    }
    tmp = str_remove(string = tmp, pattern = c("\\*"))
    tmp = str_remove(string = tmp, pattern = c("\\@"))
    tmp = str_trim(tmp)
    for (j in 1:length(tmp)){
      Plist_x$S[j,j,k] = as.numeric(tmp[j])
    }
    
    
  }
  
  # clean up Plist_x
  Plist_x$S[is.na(Plist_x$S)] = 0
  
  return(Plist_x)
}
scrutinize_tracker<-function(oneline_df){
  require(stringr)
  require(lpa.mi.src)
  
  load(with(oneline_df, paste(estwd,estfolder,estfile, sep = "/")))
  outMplus_x = list_estimates$out_Mplus
  rm(list_estimates)
  Plist_ifom = syntax2Plist(outMplus_x$input$model)
  
  Plist_est = Mplus2Qlist(params_df = outMplus_x$parameters$unstandardized)
  
  # IFOM (input from outModels) vs. est-files comparison
  L1_mu_ifom_est = max(abs(Plist_est$mu - round(Plist_ifom$mu,3)))
  L1_var_ifom_est  = max(apply(abs(Plist_est$S-round(Plist_ifom$S,3)), MARGIN=1, FUN = "max"))
  
  # .inp file vs. IFOM comparison
  # read in input model
  hi = with(oneline_df,readLines(con = paste(outwd,outfolder,paste0(outfile,".inp"), sep = "/")))
  i_0 = which(str_detect(string = hi, pattern = c("MODEL")))
  i_1 = which(str_detect(string = hi, pattern = c("OUTPUT")))
  Plist_inp = syntax2Plist(hi[seq(i_0+1,i_1-1)])
  L1_mu_inp_ifom = max(abs(Plist_inp$mu - Plist_ifom$mu))
  L1_var_inp_ifom = max(abs(Plist_inp$S - Plist_ifom$S))
  
  
  # other things to check
  error_est = length(outMplus_x$errors)
  Rcond_est = check_convergence(file = with(oneline_df, paste0(tolower(outfile),".out")), 
                                folder_wd= with(oneline_df, paste(outwd,outfolder,sep = "/")))$Rcond
  
  out_df = with(oneline_df, data.frame(nn=nn, rep=rep, z=z, pm = pm, pva=pva, m=m,kfit=kfit,data_type=data_type,
                                       converged_tracker = converged,
                                       Rcond_tracker = Rcond,
                                       L1_mu_ifom_est = L1_mu_ifom_est, 
                                       L1_var_ifom_est = L1_var_ifom_est, 
                                       L1_mu_inp_ifom = L1_mu_inp_ifom, 
                                       L1_var_inp_ifom = L1_var_inp_ifom, 
                                       error_est = error_est, 
                                       Rcond_est = Rcond_est)
  )
} 

setwd("H:/")

# Set up parallel processing 
cl<-makeSOCKcluster(Processors)
doSNOW::registerDoSNOW(cl)  


for (rep in 1:Replications){
 if (!(paste0("rep",rep,".csv") %in% list.files("S:/ping-pong"))){
   
    t1 = proc.time()
    # Serve a pong
    pong_df = data.frame(rep = rep, computer = computer_name, finish_time = NA)
    write.csv(pong_df, file = paste0("S:/ping-pong/rep",rep,".csv"), row.names = FALSE)
    
    # Copy and extract the tracker-file
    system(paste0('xcopy S:\\tracker-files\\tracker-rep',rep,'.RData H:\\tracker-files\\ /S'))
    
    # Copy over the out files
    system(paste0('xcopy S:\\out-files\\out-files-rep', rep, '.zip H:\\out-files /J /Y'))
    print("Unzipping outfile")
    system(paste0('Bandizip.exe x -y -o:H:\\out-files\\rep',rep,' H:\\out-files\\out-files-rep',rep,'.zip'))
    
    
    # Copy the est-file
    system(paste0('xcopy S:\\est-files\\rep', rep, '.zip H:\\est-files /J /Y'))
    print("Unzipping est-files")
    system(paste0('Bandizip.exe x -y -o:H:\\est-files\\rep',rep,' H:\\est-files\\rep',rep,'.zip'))
    
    # Load and clean the tracker file
    load(paste0("H:/tracker-files/tracker-rep",rep,".RData"))
    names(tracker_df)
    #table(tracker_df$outwd)
    tracker_df$outwd = "H:/out-files"
    #table(tracker_df$outfolder)
    
    #table(tracker_df$estwd)
    #table(tracker_df$estfolder)
    #table(tracker_df$estfile)
    
    # update the tracker with the new est file names
    tracker_df$estwd = "H:/est-files"
    tracker_df$estfolder = tracker_df$outfolder
    
    # Subset the tracker to the correct file
    tracker2_df = subset(tracker_df, converged == TRUE & switched==TRUE)
    
    pb <- pbapply::timerProgressBar(max = nrow(tracker2_df), style = 1, width = getOption("width")/4)
    progress <- function(x){setTimerProgressBar(pb, x)}
    opts <- list(progress = progress)
    tmplist<-foreach(x = 1:nrow(tracker2_df),
                     .packages = c("lpa.mi.src","stringr"),
                     .inorder = TRUE, .options.snow = opts) %dopar% {
                       
                       outline_x = scrutinize_tracker(tracker2_df[x,])
                       return(outline_x)
                       
                     }
    
    # Create and save the scrutiny file
    scrutinize_df = data.table::rbindlist(tmplist)
    save(scrutinize_df, file = paste0("H:/scrutiny-files/scrutiny-table-rep", rep,".RData"))
    system(paste0('xcopy H:\\scrutiny-files\\scrutiny-table-rep',rep,'.RData S:\\scrutiny-files /Y /J'))
    
    # Clean up the H drive
    system('rm -r H:\\tracker-files')
    system('mkdir H:\\tracker-files')
    system('rm -r H:\\est-files')
    system('mkdir H:\\est-files')
    system('rm -r H:\\out-files')
    system('mkdir H:\\out-files')
    
    # Update pong with finish time
    toc = proc.time()-t1
    pong_df$finish_time = toc[3]
    write.csv(pong_df, file = paste0("S:/ping-pong/rep",rep,".csv"), row.names = FALSE)

 }
}

close(pb)
stopCluster(cl)

