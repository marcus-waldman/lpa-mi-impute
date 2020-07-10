rm(list = ls())
gc()

library(plyr)
library(tidyverse)
library(data.table)
library(pbapply)
library(stringr)
library(lpa.mi.src)
library(tangram)
library(dplyr)
library(qwraps2)
library(Hmisc)
library(flexclust)

library(doParallel)
library(foreach)
library(doRNG)
require(snow)
require(doSNOW)
require(foreach)
require(pbapply)



# Directories
dropbox_wd = "D:/Dropbox"
#dropbox_wd = "C:/Users/marcu/Dropbox"
results_wd = paste0(dropbox_wd, "/Dissertation/lpa-mi-impute/stage4c-combine-results")
stage6_wd = paste0(dropbox_wd, "/Dissertation/lpa-mi-impute/stage6-classification-accuracy")
environment_wd = paste0(dropbox_wd,"/Dissertation/environmental-variables/")


setwd(environment_wd)
load(file = "labels_for_conditions.Rdata")
load(file = "data-conditions.Rdata")


# # Get classification
# # Load in the data
# setwd(results_wd)
# load(file ="summaries-combined-results-lpa-mi-impute.RData")
# 
# # Append in the classification data
# system('rm -r H:\\classification-accuracy')
# system('mkdir H:\\classification-accuracy')
# system('xcopy S:\\classify-accuracy-files\\* H:\\classification-accuracy\\ /J /Y')
# files_vec = list.files("H:/classification-accuracy/", full.names = T)
# classify_df<-
#   lapply(X = 1:500, FUN = function(x){load(files_vec[x]); return(out_x %>% transform(rep = x))}) %>% 
#   rbindlist
# classify_df = classify_df %>% 
#   transform(Imputation_Procedure = mapvalues(data_type, 
#                                              from = c("Complete data", "Observed data", "Imputation"), 
#                                              to = c("No Missingness", "FIML", "NA"))
#   )
# classify_df$Imputation_Procedure[classify_df$pva==1] = "Amelia"
# classify_df$Imputation_Procedure[classify_df$pva==2] = "PMM"
# classify_df$Imputation_Procedure[classify_df$pva==3] = "CART"
# classify_df$Imputation_Procedure[classify_df$pva==4] = "RF"
# 
# 
# 
# 
# # Merge with summary
# summaries_combined_df =  summaries_combined_df %>% filter(Imputation_Procedure != "Pattern Mixture")
# nrow(summaries_combined_df)
# summaries_combined_df = summaries_combined_df %>% left_join(classify_df %>% select(-c("data_type")) %>% data.frame(), by = c("rep","z","pva","pm","Imputation_Procedure"))
# nrow(summaries_combined_df)
# 
# # Processors = 10
# # cl<-makeSOCKcluster(Processors)
# # doSNOW::registerDoSNOW(cl)   
# # 
# # 
# # pb <- pbapply::timerProgressBar(max = nrow(summaries_combined_df), style = 1, width = getOption("width")/4)
# # progress <- function(x){setTimerProgressBar(pb, x)}
# # opts <- list(progress = progress)
# # tmp_list<-
# #   foreach(x = 1:nrow(summaries_combined_df),
# #           .packages = c("tidyverse","flexclust"), 
# #           .inorder = TRUE, 
# #           .options.snow = opts) %dopar% {
# # 
# # 
# # 
# #                hi = summaries_combined_df[x, ] %>% select(starts_with("kappa")) %>% as.integer() %>% matrix(nrow = 3)
# #                if(sum(is.na(hi))==0){
# #                 
# #                  Rand = randIndex(hi %>% as.table(), correct = FALSE)
# # 
# #                } else {
# #                  Rand = NA
# #                }
# #                return(with(summaries_combined_df, 
# #                           data.frame(Imputation_Procedure = Imputation_Procedure[x],
# #                                      rep = rep[x],
# #                                      z = z[x],
# #                                      pm = pm[x], 
# #                                      pva = pva[x], 
# #                                      Rand = Rand)))
# #           }
# # 
# # 
# # 
# #
# tmp_df = rbindlist(tmp_list)
# summaries_v2_combined<- summaries_combined_df %>% left_join(tmp_df %>% select(c("rep","z","pva","pm","Imputation_Procedure","Rand")), 
#                                    by = c("rep","z","pva","pm","Imputation_Procedure"), 
#                                    )
# setwd(results_wd)
# save(summaries_v2_combined, file = "summaries-v2-combined-results-rand-lpa-mi-impute.RData")

# stopCluster(cl)

setwd(results_wd)
load("summaries-v2-combined-results-rand-lpa-mi-impute.RData")

agg_df = summaries_v2_combined %>% 
         left_join(conditions_labels, by = "z") %>% 
         group_by(Imputation_Procedure,Mixing,Sample.Size,Separation) %>% 
         summarise(Rand = mean(Rand, na.rm = T)) %>% 
         filter(Imputation_Procedure != "PMM")
ggplot(agg_df %>% filter(startsWith(as.character(Sample.Size), "Small")), aes(x = Imputation_Procedure, y = Rand)) + geom_point() + facet_grid(Separation~Mixing)

ggplot(agg_df %>% filter(startsWith(as.character(Sample.Size), "Large")), aes(x = Imputation_Procedure, y = Rand)) + geom_point() + facet_grid(Separation~Mixing)


# Long story short: The improvement in the Rand Index is not impressive, but non the less performs better than FIML alwayhs and better than Amelia.