rm(list = ls())

library(plyr)
library(tidyverse)
library(data.table)
library(pbapply)


# Copy over the results


results_wd = "S:/results"
stage_wd = "D:/Dropbox/Dissertation/lpa-mi-impute/stage5-analyze-results"

load("S:/environment-mi-impute-stage0 Jan 06 2019 13 29.RData")


all_vec = list.files(path = results_wd)
parameter_files = all_vec[startsWith(all_vec, "results-parameters")]
summaries_files = all_vec[startsWith(all_vec, "results-summaries")]



# # Look at parameters
tmp<-pblapply(X = seq(1,length(parameter_files)),
              FUN = function(x){load(file = paste0(results_wd,"/",parameter_files[x]));return(parameters_df)})
parameters_df = rbindlist(tmp)


# Is there bias in FIML if no correlation, 

# Means
means_df = subset(parameters_df,
               paramHeader=="Means" & startsWith(param,"Y") & Imputation_Procedure != "No Missingness")


table(subset(parameters_df, Sample_Size == "Small N" & Association == "No Correlation")$Separation)
table(subset(parameters_df, Sample_Size == "Small N" & Association == "No Correlation")$Mixing)
table(subset(parameters_df, Sample_Size == "Small N" & Association == "No Correlation")$Mean_Differences)
table(subset(parameters_df, Sample_Size == "Small N" & Association == "No Correlation")$Moderation)


ggplot(subset(means_df,Sample_Size == "Large N"), aes(deviance)) + geom_density() + 
  facet_grid(Imputation_Procedure~LatentClass) + 
  coord_cartesian(xlim = c(-0.5,0.5))

pdf(file = paste0(stage_wd,"/FIML_LargeN.pdf"), width = 8.5, height = 22)
ggplot(subset(means_df,Sample_Size == "Large N" & Imputation_Procedure == "FIML"), aes(deviance)) + geom_density() + 
  facet_grid(Mixing+Moderation+Separation+Mean_Differences~LatentClass) + 
  coord_cartesian(xlim = c(-0.5,0.5)) + 
  ggtitle(label = "FIML, Large N")
dev.off()

pdf(file = paste0(stage_wd,"/Amelia_LargeN.pdf"), width = 8.5, height = 22)
ggplot(subset(means_df,Sample_Size == "Large N" & Imputation_Procedure == "Amelia"), aes(deviance)) + geom_density() + 
  facet_grid(Mixing+Moderation+Separation+Mean_Differences~LatentClass) + 
  coord_cartesian(xlim = c(-0.5,0.5)) + 
  ggtitle(label = "Amelia, Large N")
dev.off()

# Explore FIML
ggplot(subset(means_df,Imputation_Procedure == "FIML" & Sample_Size == "Large N"), aes(deviance)) + 
  geom_density() + facet_grid(Association~LatentClass) + 
  coord_cartesian(xlim = c(-0.5,0.5))

ggplot(subset(means_df,Sample_Size == "Small N"), aes(deviance)) + geom_density() + facet_grid(Imputation_Procedure~LatentClass) + 
  coord_cartesian(xlim = c(-0.5,0.5))



# One-Way-Anova: (a) Separation, (b) Mean Differences, (c) Association, (d) Mixing, (e) Moderation
bias_means_ssize = ddply(means_df, .(Sample_Size,Imputation_Procedure,LatentClass), summarize,
                         bias = mean(deviance), 
                         N = length(deviance),
                         sd_bias = sd(deviance))
ggplot(bias_means_ssize, aes(x = LatentClass, y = bias)) + 
  geom_text(aes(label = Imputation_Procedure)) + 
  facet_grid(Sample_Size~.) + 
  geom_hline(yintercept = 0) 


# Two-Way Anova plot across: (a) Separation, (b) Mean Differences, (c) Association, (d) Mixing, (e) Moderation
bias_means_ssize2 = ddply(means_df, .(Sample_Size,Imputation_Procedure,LatentClass, Mixing, Separation), summarize,
                          bias = mean(deviance), 
                          N = length(deviance),
                          sd_bias = sd(deviance))
ggplot(subset(bias_means_ssize2, Sample_Size == "Large N"), aes(x = LatentClass, y = bias)) + 
  geom_text(aes(label = Imputation_Procedure)) + 
  facet_grid(Sample_Size~.) + 
  geom_hline(yintercept = 0) + 
  facet_grid(Mixing~Separation)


bias_means_ssize2 = ddply(means_df, .(Sample_Size,Imputation_Procedure,LatentClass, Mixing, Separation), summarize,
                          bias = mean(deviance), 
                          N = length(deviance),
                          sd_bias = sd(deviance))
ggplot(subset(bias_means_ssize2, Sample_Size == "Large N"), aes(x = LatentClass, y = bias)) + 
  geom_text(aes(label = Imputation_Procedure)) + 
  facet_grid(Sample_Size~.) + 
  geom_hline(yintercept = 0) + 
  facet_grid(Mixing~Separation)




ggplot(subset(bias_means_ssize2, Sample_Size == "Large N"), aes(x = LatentClass, y = bias)) + 
  geom_text(aes(label = Imputation_Procedure)) + 
  facet_grid(Sample_Size~.) + 
  geom_hline(yintercept = 0) + 
  facet_grid(Mixing~Separation)


ggplot(subset(bias_means_ssize, Sample_Size == "Small N"), aes(x = LatentClass, y = bias)) + 
  geom_text(aes(label = Imputation_Procedure)) + 
  facet_grid(Sample_Size~.) + 
  geom_hline(yintercept = 0) 

ggplot(subset(bias_means_ssize2, Sample_Size == "Small N"), aes(x = LatentClass, y = bias)) + 
  geom_text(aes(label = Imputation_Procedure)) + 
  facet_grid(Sample_Size~.) + 
  geom_hline(yintercept = 0) + 
  facet_grid(Mixing~Separation)

  


# Summaries

tmp<-pblapply(X = seq(1,length(parameter_files)),
              FUN = function(x){load(file = paste0(results_wd,"/",summaries_files[x]));return(summaries_df)})
summaries_df = rbindlist(tmp)

KL_df = ddply(summaries_df, .(Sample_Size,Imputation_Procedure), summarize,
              KL = weighted.mean(KL, se_KL^(-2)), 
              N = length(KL),
              sd_KL = sd(KL))

ggplot(subset(KL_df, Sample_Size == "Large N"), aes(x = Imputation_Procedure, y = KL)) + 
  geom_point() + 
  facet_grid(Sample_Size~.) 


ggplot(subset(KL_df, Sample_Size == "Small N"), aes(x = Imputation_Procedure, y = KL)) + 
  geom_point() + 
  facet_grid(Sample_Size~.) 




MD_df = rbindlist(tmp_list)


pick1_df = left_join(x=pick1_df,y=MD_df)
ddply(pick1_df, .(Imputation_Procedure), summarize, 
      MD = mean(bias_MD))
