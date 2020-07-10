rm(list = ls())

library(plyr)
library(tidyverse)
library(data.table)
library(pbapply)


# Copy over the results


results_wd = "S:/results"
stage_wd = "D:/Dropbox/Dissertation/lpa-mi-impute/stage5-analyze-results"



all_vec = list.files(path = results_wd)
parameter_files = all_vec[startsWith(all_vec, "results-parameters")]
summaries_files = all_vec[startsWith(all_vec, "results-summaries")]


# # Look at parameters
tmp<-pblapply(X = seq(1,length(parameter_files)),
              FUN = function(x){load(file = paste0(results_wd,"/",parameter_files[x]));return(parameters_df)})
parameters_df = rbindlist(tmp)
# rm(tmp)
# 
# # 
# # # Means
# # 
# # tmp_df = subset(parameters_df, 
# #                paramHeader=="Means" & startsWith(param,"Y") & Sample_Size == "Large N" & Imputation_Procedure != "CART") %>% 
# #          ddply(.(z,LatentClass,param,Imputation_Procedure),summarize, bias = mean(deviance))
# # 
# # cart_df =  subset(parameters_df, 
# #                   paramHeader=="Means" & startsWith(param,"Y") & Sample_Size == "Large N" & Imputation_Procedure == "CART") %>% 
# #   ddply(.(z,LatentClass,param,Imputation_Procedure),summarize, bias = mean(deviance))
# # 
# # ggplot() + geom_point(data = tmp_df, aes(x = z, y = bias, col = Imputation_Procedure)) + 
# #            geom_point(data = cart_df, aes(x = z, y = bias)) + 
# #            facet_grid(param~LatentClass, scales = "free")
# #            
# # 
# # # Variances
# # 
# # tmp_df = subset(parameters_df, 
# #                 paramHeader=="Variances" & startsWith(param,"Y") & Sample_Size == "Large N" & Imputation_Procedure != "CART") %>% 
# #   ddply(.(z,LatentClass,param,Imputation_Procedure),summarize, bias = mean(deviance))
# # 
# # cart_df =  subset(parameters_df, 
# #                   paramHeader=="Variances" & startsWith(param,"Y") & Sample_Size == "Large N" & Imputation_Procedure == "CART") %>% 
# #   ddply(.(z,LatentClass,param,Imputation_Procedure),summarize, bias = mean(deviance))
# # 
# # ggplot() + geom_point(data = tmp_df, aes(x = z, y = bias, col = Imputation_Procedure)) + 
# #   geom_point(data = cart_df, aes(x = z, y = bias)) + 
# #   facet_grid(param~LatentClass, scales = "free")
# # 
# # 
# # 
# # # Class proportions
# # tmp_df = subset(parameters_df, 
# #                 paramHeader=="Means" & startsWith(param,"C#") & Sample_Size == "Large N" & Imputation_Procedure != "CART") %>% 
# #   ddply(.(z,param,Imputation_Procedure,Mixing),summarize, bias = mean(deviance))
# # 
# # cart_df =  subset(parameters_df, 
# #                   paramHeader=="Means" & startsWith(param,"C#") & Sample_Size == "Large N" & Imputation_Procedure == "CART") %>% 
# #   ddply(.(z,param,Imputation_Procedure,Mixing),summarize, bias = mean(deviance))
# # 
# # ggplot() + geom_point(data = tmp_df, aes(x = z, y = abs(bias), col = Imputation_Procedure)) + 
# #   geom_point(data = cart_df, aes(x = z, y = abs(bias))) + 
# #   facet_grid(param~Mixing, scales = "free")
# # 
# 
# 
# # Summaries
# 
# 
# 
tmp<-pblapply(X = seq(1,length(parameter_files)),
              FUN = function(x){load(file = paste0(results_wd,"/",summaries_files[x]));return(summaries_df)})
summaries_df = rbindlist(tmp)

#  tmp_df = subset(summaries_df, Imputation_Procedure != "CART") %>% 
#           ddply(.(z,Imputation_Procedure),summarize, KL = weighted.mean(KL,se_KL^(-2)))
#  
#  CART_df = subset(summaries_df, Imputation_Procedure == "CART") %>% 
#    ddply(.(z,Imputation_Procedure),summarize, KL = weighted.mean(KL,se_KL^(-2)))
#  
# 
#   ggplot() + geom_point(data = tmp_df, aes(x = z, y = KL, col = Imputation_Procedure)) + 
#     geom_point(data = CART_df, aes(x = z, y = KL))
#   
# 
# summaries_df$Imputation_Procedure = as.factor(summaries_df$Imputation_Procedure)
# summaries_df$Imputation_Procedure = relevel(summaries_df$Imputation_Procedure, ref = "No Missingness")
# 
# obj_lm = lm(KL~Imputation_Procedure, data = summaries_df, weight = (summaries_df$se)^(-2))
# 


# Class separation
library(lpa.mi.src)
pick1_df = parameters_df %>% 
           filter(paramHeader=="Means" & param=="Y1" & LatentClass==1) %>% 
           transform(id = 1:nrow(.)) %>% 
           select(rep,z,pm,kfit,Imputation_Procedure,id)
nrow(parameters_df)
parameters_df = parameters_df %>% left_join(y = pick1_df)
nrow(parameters_df)

max(pick1_df$id)

load("S:/environment-mi-impute-stage0 Jan 06 2019 13 29.RData")

library(parallel)
cl<-makeCluster(10)
clusterExport(cl, c("parameters_df", "data_conditions"))


# library(pbapply)
# tmp_list = pblapply(X = 1:nrow(pick1_df), 
#                     FUN = function(x){
#                           require(lpa.mi.src)
#                           require(stats)
#                           parameters_x = subset(parameters_df, id == x)
#                           z_x = parameters_x$z[1]
#                           Qlist_x = Mplus2Qlist(parameters_x)
#                           J_Y_x = data_conditions$J_Y[z_x]
#                           K_x = data_conditions$K[z_x]
#                           MD_mat = mat.or.vec(nr = K_x, nc = K_x) + NA
#                           for (i in seq(2,K_x)){
#                             for (j in seq(1,i-1)){
#                               MD_mat[i,j]<-MD_mat[j,i] <- with(Qlist_x, stats::mahalanobis(x = mu[,i], center = mu[,j], cov = diag(J_Y_x)))
#                             }
#                           }
#                           MD_mean  = mean(as.vector(MD_mat), na.rm = TRUE)
#                           return(data.frame(id = x, bias_MD = MD_mean - data_conditions$MD[z_x]))
#                     }, 
#                     cl = cl)
# stopCluster(cl)
# save(tmp_list, file = paste0(results_wd,"/MH_list.RData"))

MD_df = rbindlist(tmp_list)


pick1_df = left_join(x=pick1_df,y=MD_df)
ddply(pick1_df, .(Imputation_Procedure), summarize, 
      MD = mean(bias_MD))
