rm(list = ls())
library(Hmisc)
library(tidyverse)
library(dplyr)


load("D:/Dropbox/Dissertation/environmental-variables/environment-mi-impute-stage0 Jan 06 2019 13 29.Rdata")
rm(list = as.character(setdiff(ls(),"data_conditions")))

data_conditions =  transform(data_conditions,
                             z = 1:nrow(data_conditions),
                             Mixing = ordered(ifelse(class_size == "A",0, 1), levels = c(0,1), label = c("Equal Mix.", "Unequal Mix.")),
                             Sample.Size = ordered(ifelse(N==300, 0,1), levels = c(0,1), label = c("Small Sample (N=300)", "Large Sample (N=1200)")), 
                             Separation = ordered(ifelse(MD==2.87,0,1), levels = c(0,1), label = c("Weakly Separated", "Strongly Separated")),  
                             AV.Correlation = ordered(ifelse(rho_YX==0.4,0,1), levels = c(0,1), label = c("AV-Strong Corr.", "AV-No Corr.")), 
                             AV.Mean.Diff = ordered(ifelse(dX == TRUE,0,1),levels = c(0,1), label = c("AV-Mean Class Diff.", "AV-No Mean Class Diff.")), 
                             AV.Moderation = ordered(ifelse(C_modifies_YX==TRUE,0,1), levels = c(0,1) ,c("AV-Moderation", "AV-No Moderation")))
conditions_labels = data_conditions %>% 
                  select(z,Mixing, Sample.Size,Separation,AV.Correlation,AV.Mean.Diff,AV.Moderation)

setwd("D:/Dropbox/Dissertation/environmental-variables/")

save(conditions_labels, file = "labels_for_conditions.Rdata")
