#  .rs.restartR()

rm(list = ls())

library(plyr)
library(tidyverse)
library(lpa.mi.src)
library(doParallel)
library(foreach)
library(doRNG)


Replications = 500
RndSeed = 42
this_file = "mi-impute-stage0-z1to40-Gaby.R"

# Acer-laptop
  #dropbox_wd =  "C:/Users/marcu/Dropbox"  #  Dropbox folder on local machine
  #datfiles_dir = "C:/Users/marcu/Documents/lpa-mi-create-stage0" #Stores Mplus-read dat files
  #Rdata_dir = "C:/Users/marcu/Documents/lpa-mi-create-stage0"
# Queen Mary's Revenge
  #dropbox_wd =  "D:/Dropbox"  #  Dropbox folder on local machine
  #datfiles_dir = #Stores Mplus-read dat files
# Xi-GSU
  #dropbox_wd = "C:/Users/mwaldman1/Dropbox"
  #datfiles_dir =  #Stores Mplus-read dat files
# Big Bertha
  #dropbox_wd =  "D:/Dropbox/Dropbox"  #  Dropbox folder on local machine
  #datfiles_dir = #Stores Mplus-read dat files
  #Rdata_dir = 
# Gaby
  dropbox_wd = "D:/Dropbox"
  datfiles_dir = "E:" #Stores Mplus-read dat files
  Rdata_dir = "E:"

# Fully-crossed simulation conditions
pctmiss_vec = c(0.5); pm_vec = c(1)  #  Percent of observations with *at least* one missing value (indicator or missing data correlate) -> if other specification needed one must modify list_inds.miss generation in get_obs_data.R
N_vec = c(300, 1200)          #  Sample sizes within each replication (small - 300, large - 1200)
J_Y_vec = c(4)                #  Number of latent class indicators (K = 4)
J_Xcom_vec = c(1)             #  Number of complete data missing data correlates
J_Xinc_vec = c(0)             #  Number of missing data correlates which themselves contain missing data
MD_vec = c(2.87,3.7)          #  Mahalanobis disance separating classes (low - 2.87, high = 3.70)
pi_list = list(               #  Mixing and (by proxy) the number of classes, K (maximum of 3)
  A = rep(1/3,3), 
  B = c(0.45,0.45,0.1) 
)
rho_YX_vec = c(0.4, 0)          #  Strength of missing data correlates, correlation
dX_vec = c(TRUE, FALSE)         #  Between class differences in missing data correlates, Cohen's d
C_modifies_YX_vec = c(FALSE, TRUE)    #  Whether or not the Y-X relationships are different between latent classes 
t_rotate_vec = c(0)             #  Angle for rotation matrix for LPA indicators 


######################################################################################################################
#      Pre-process the data generating conditions                                                                    #
######################################################################################################################
# Identify all conditions to be evaluated across reps = 1:Replications for each method in method_list 
data_conditions = data.frame(expand.grid(N = N_vec, J_Y = J_Y_vec, J_Xcom = J_Xcom_vec, 
                                         J_Xinc = J_Xinc_vec, MD = MD_vec, class_size = names(pi_list), 
                                         rho_YX = rho_YX_vec, dX = dX_vec,
                                         C_modifies_YX = C_modifies_YX_vec, t_rotate = t_rotate_vec))
# Remove the AV as non-confounder conditions
z_remove = with(data_conditions, which(rho_YX==0 & dX == FALSE | rho_YX==0 & C_modifies_YX_vec==TRUE))
if(length(z_remove)>0){data_conditions = data_conditions[-z_remove, ]}
# Add in the number of classes corresponding to each condition and the marginal probabilities (pi_k)
Kmax = max(simplify2array(lapply(pi_list, FUN = length)))
pi_df = data.frame(mat.or.vec(nr = length(pi_list), nc = Kmax + 2))+NA
names(pi_df) = c("class_size","K", paste("pi_", 1:Kmax, sep = ""))
for(p in 1:nrow(pi_df)){
  pi_df$class_size[p] = names(pi_list)[p]
  pi_df$K[p] = length(pi_list[[p]])
  pi_df[p,seq(3,pi_df$K[p]+2)] = pi_list[[p]]
}
data_conditions = merge(x = data_conditions, y = pi_df)
rm(pi_df)
# Ensure that there is at least one data condition
Ztot = nrow(data_conditions)
if (Ztot==0){stop("Total number supported of data conditions is zero.")}
#Sort and check pctmiss is between 0 and 1
pctmiss_vec = sort(pctmiss_vec)
if( (sum(pctmiss_vec<=0) + sum(pctmiss_vec>=1))>0){stop("elements in pctmiss_vec must be between 0 and 1 (exclusive).")}
# Define working directory
main_wd = paste(dropbox_wd, "/Dissertation/lpa-mi-impute", sep = "")
# Initiate a list of folder files
dffolderfiles = NULL



######################################################################################################################
#      Create folders and subfolders in temporary directory for read-ins and write-outs at each replication          #
######################################################################################################################
# Create folder name for temporary folder (saved in results and in temp files)
temp_fname = paste("mi-impute-stage0", format(Sys.time(), "%b %d %Y %H %M"), sep = " ")
# Create temporary folder
temp_wd =  paste(datfiles_dir, temp_fname, sep = "/")
dir.create(temp_wd)
setwd(temp_wd)
# Create replication-specific temporary folder
datfiles_wd_rep_vec = NULL
for (rep in 1:Replications){
  temp_wd_rep = paste0(temp_wd,"/rep",rep, sep = "")
  dir.create(temp_wd_rep)
  datfiles_wd_rep_vec = c(datfiles_wd_rep_vec, temp_wd_rep)
}
# Within the replication-specific temporary folder, create subfolders by z-th conditions
for (rep in 1:Replications){
    for(z in 1:Ztot){
    dir.create(paste0(datfiles_wd_rep_vec[rep],"/z",z))
    dir.create(paste0(datfiles_wd_rep_vec[rep],"/z",z,"/Complete data"))
    dir.create(paste0(datfiles_wd_rep_vec[rep],"/z",z,"/Observed data"))
    for (i in 1:length(pctmiss_vec)){
      dir.create(paste0(datfiles_wd_rep_vec[rep],"/z",z,"/Observed data/pm", pm_vec[i]))
    }
  }
}



######################################################################################################################
#    Create file that stores results at each replication eventually to be zipped and copied to dropbox results       #
######################################################################################################################
lists_wd = paste0(Rdata_dir,"/Rdata-",temp_fname)
dir.create(lists_wd)
# Save environment
save.image(file = paste0(lists_wd,"/environment-", temp_fname,".RData"))
# Copy all R files in main directory
file.copy(from = paste0(main_wd,"/",this_file), to = paste0(lists_wd,"/","executed-R-code ", temp_fname, ".R"))




######################################################################################################################
#      Initiate other variables and write out the environment                                                        #
######################################################################################################################
set.seed(RndSeed)
seeds_vec = .Random.seed[1:Replications]
pb <- winProgressBar(title = "progress bar", min = 0,
                     max = Replications, width = 300)
t1 = proc.time()      
#rep = 1
#z = 1
sink(file = paste0(lists_wd, "/sink-stage0 rep", rep,".txt"))
for(rep in 1:Replications){
  RndSeed_rep = seeds_vec[rep]
  set.seed(RndSeed_rep)
  for(z in 1:Ztot){
    print(paste0("--------------- rep = ", rep, ",  z = ", z, "  ---------------"))
    
    pop_params_z = lpa.mi.src::get_FMM_params(z = z, data_conditions = data_conditions)
    
    
    # Complete Data
    list_complete = lpa.mi.src::get_complete_data(z = z, 
                                                  data_conditions = data_conditions,
                                                  rep = rep, 
                                                  save_it = TRUE, 
                                                  temp_wd_rep_vec = datfiles_wd_rep_vec)
    print("  >> list_complete generated")
    
    # Observed data
    list_observed = lpa.mi.src::get_obs_data(z=z, 
                                             df =  list_complete$dfcom,
                                             pctmiss_vec = pctmiss_vec,
                                             data_conditions = data_conditions,
                                             rep = rep,
                                             save_it = TRUE, 
                                             temp_wd_rep_vec = datfiles_wd_rep_vec,
                                             pm_vec = pm_vec)
    print("  >> list_observed generated")
    
    # Details descriptves
    print("  >> observed data rates")
    print(list_observed$obs_rates)
    print("  >> observed data counts")
    print(list_observed$obs_counts)
    
    save(list = c("list_complete","list_observed"), 
         file = paste0(lists_wd, "/list-stage0 rep", rep, " z",z,".RData"))
    print("  >> RData file saved")
    print(" ")
    print(" ")
    
    
  } ## end for z = 1:Ztot
  setWinProgressBar(pb, rep, title=paste( round(rep/Replications*100, 0),"% done"))
} #for rep in 1:Replications
sink()

proc.time()-t1   
close(pb)




