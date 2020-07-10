rm(list = ls())

require(pbapply)
require(plyr)
require(tidyverse)
require(data.table)
library(psych)
library(gmodels)


# Load in the scrutiny file
files_scrutiny = list.files("S:/scrutiny-files/", full.names = TRUE)


list_scrutiny<-pblapply(X = files_scrutiny, 
                        FUN = function(x){load(x); return(scrutinize_df)})
scrutiny_df = rbindlist(list_scrutiny)
rm(list_scrutiny)
gc()


# Load in the tracker files
files_tracker = list.files("S:/tracker-files/", full.names = TRUE, pattern = ".RData")
list_trackers = pblapply(X = files_tracker, 
                         FUN = function(x){load(x);return(tracker_df)})
tracker_df = rbindlist(list_trackers)
rm(list_trackers)
gc()


# Summarize the scrutiny file
head(scrutiny_df)
scrutiny_df =  scrutiny_df %>% 
  transform(problem_inp_ifom = (L1_mu_inp_ifom>=1E-5 | L1_var_inp_ifom>=1E-5),
            problem_ifom_est = (L1_mu_ifom_est>=1E-5 | L1_var_ifom_est>=1E-5),
            problem_Rcond = Rcond_est<1E-6)

#with(scrutiny_df, gmodels::CrossTable(x = problem_inp_ifom, y = problem_ifom_est, prop.r = T, prop.c = T, prop.t = T, prop.chisq = F))
#with(scrutiny_df, gmodels::CrossTable(x = problem_Rcond, y = problem_ifom_est, prop.r = T, prop.c = T, prop.t = T, prop.chisq = F))
#with(scrutiny_df, gmodels::CrossTable(x = problem_Rcond, y = problem_inp_ifom, prop.r = T, prop.c = T, prop.t = T, prop.chisq = F))


# Define all that didn't pass scrutiny
scrutiny_df$pass_scrutiny = TRUE
table(scrutiny_df$problem_inp_ifom, useNA = "always")
scrutiny_df$pass_scrutiny[which(scrutiny_df$problem_inp_ifom)] = FALSE

table(scrutiny_df$problem_ifom_est, useNA = "always")
scrutiny_df$pass_scrutiny[which(scrutiny_df$problem_ifom_est)] = FALSE

table(scrutiny_df$Rcond_est<1E-6, useNA = "always")
scrutiny_df$pass_scrutiny[which(scrutiny_df$Rcond_est<1E-6)] = FALSE
sum(is.na(scrutiny_df$Rcond_est))
scrutiny_df$pass_scrutiny[which(is.na(scrutiny_df$Rcond_est))] = FALSE


table(scrutiny_df$error_est, useNA = "always")
scrutiny_df$pass_scrutiny[which(scrutiny_df$error_est>0)]  = FALSE


# Merge in pass scrutiny
nrow(tracker_df)
tracker_df = tracker_df %>% left_join(scrutiny_df[,c("nn","rep","pass_scrutiny")]) 
nrow(tracker_df)

# common sense quick checks
nrow(scrutiny_df) 
sum(!is.na(tracker_df$pass_scrutiny)) #looks good

sum(scrutiny_df$pass_scrutiny)
sum(tracker_df$pass_scrutiny, na.rm = TRUE) # PASS

with(tracker_df, table(switched, pass_scrutiny, useNA = "always")) #looks good

# Now save the individual trackers to 

table(scrutiny_df$pass_scrutiny, useNA = "always")

rm("failed_df")
rm("scrutiny_df")
rm(files_scrutiny)
rm(files_tracker)
gc()

tracker_all = tracker_df
rm(tracker_df)

load(file = "S:/tracker-files/tracker-rep1.RData")

for (j in 1:32){
  print(names(tracker_df)[j])
  print(class(tracker_df[,j]))
  print(class(tracker_all[,j]))
  print(" ")
  print(" ")
  
}

for (j in 1:32){
  if (is.factor(tracker_all[,j])){
    print(names(tracker_all)[j])
    tracker_all[,j] = as.character(tracker_all[,j])
  }
}


tracker_all$rep = as.integer(tracker_all$rep)

rm(tracker_df)

library(parallel)
list_save<-pblapply(X = 1:500, 
                    FUN = function(x){
                      tracker_df = subset(tracker_all, rep == x);
                      save(tracker_df, file = paste0("S:/tracker-files/tracker-rep",x,".RData"), compress = TRUE);
                      return(NULL)
                    }, 
                    cl = NULL)

load('S:/tracker-files/archieve/tracker-rep1.RData')
tracker1_orig = tracker_df
rm(tracker_df)

load('S:/tracker-files/tracker-rep1.RData')
tracker1_new = tracker_df
rm(tracker_df)

