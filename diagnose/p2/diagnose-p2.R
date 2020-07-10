rm(list = ls())

setwd("C:/Users/marcu/Dropbox/Dissertation/lpa-mi-impute/diagnose/p2")

load("dff_imputed.RData")
load("environment.RData")

library(lpa.mi.src)
library(plyr)
library(tidyverse)
library(lpa.mi.src)
library(doParallel)
library(foreach)
library(doRNG)
library(MplusAutomation)


p = 2
z = 1
pm = 1
pva = 1

temp_wd_p_vec = c(getwd(), getwd())

pop_params_z = get_FMM_params(z, data_conditions)

out = fit_and_do_rubinrules(dff_imputed, methods_list, pop_params_z,
                                data_conditions, temp_wd_p_vec, p, z, pm, pva)