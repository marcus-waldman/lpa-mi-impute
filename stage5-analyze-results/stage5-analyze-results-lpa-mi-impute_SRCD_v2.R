rm(list = ls())

library(plyr)
library(tidyverse)
library(data.table)
library(pbapply)
library(stringr)


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
                  paramHeader=="Means" & startsWith(param,"Y"))

pdf(file = paste0(stage_wd,"/Means-Density-LargeN-Correlation.pdf"), width = 8.5, height = 22)
methods_vec = unique(means_df$Imputation_Procedure)
for(i in 1:length(methods_vec)){
  df = subset(means_df, Imputation_Procedure == as.character(methods_vec[i]) & Sample_Size == "Large N" & Association == "Correlation")
  aggdf = ddply(df,.(Mixing,Moderation,Separation,Mean_Differences,LatentClass), summarize, bias = mean(deviance, na.rm = TRUE))
  p = ggplot(df, aes(deviance)) + 
    geom_histogram(bins=50, fill = "blue", alpha = 0.5) +
    geom_vline(xintercept = 0, col = "black", size = 1) + 
    geom_vline(data = aggdf, aes(xintercept = bias), col = "blue", linetype = "dashed") +
    facet_grid(Mixing+Moderation+Separation+Mean_Differences~LatentClass) +
    coord_cartesian(xlim = c(-0.5,0.5)) + 
    ggtitle(label = paste0(methods_vec[i],", Large N with Correlation")) 
  print(p)
}
dev.off()

pdf(file = paste0(stage_wd,"/Means-Density-SmallN-Correlation.pdf"), width = 8.5, height = 22)
methods_vec = unique(means_df$Imputation_Procedure)
for(i in 1:length(methods_vec)){
  df = subset(means_df, Imputation_Procedure == as.character(methods_vec[i]) & Sample_Size == "Small N" & Association == "Correlation")
  aggdf = ddply(df,.(Mixing,Moderation,Separation,Mean_Differences,LatentClass), summarize, bias = mean(deviance, na.rm = TRUE))
  p = ggplot(df, aes(deviance)) + 
    geom_histogram(bins=50, fill = "blue", alpha = 0.5) +
    geom_vline(xintercept = 0, col = "black", size = 1) + 
    geom_vline(data = aggdf, aes(xintercept = bias), col = "blue", linetype = "dashed") +
    facet_grid(Mixing+Moderation+Separation+Mean_Differences~LatentClass) +
    coord_cartesian(xlim = c(-1,1)) + 
    ggtitle(label = paste0(methods_vec[i],", Small N with Correlation")) 
  print(p)
}
dev.off()

pdf(file = paste0(stage_wd,"/Means-Density-SmallN-NoCorrelation.pdf"), width = 8.5, height = 11)
methods_vec = unique(means_df$Imputation_Procedure)
for(i in 1:length(methods_vec)){
  df = subset(means_df, Imputation_Procedure == as.character(methods_vec[i]) & Sample_Size == "Small N" & Association == "No Correlation")
  aggdf = ddply(df,.(Mixing,Moderation,Separation,Mean_Differences,LatentClass), summarize, bias = mean(deviance, na.rm = TRUE))
  p = ggplot(df, aes(deviance)) + 
    geom_histogram(bins=50, fill = "blue", alpha = 0.5) +
    geom_vline(xintercept = 0, col = "black", size = 1) + 
    geom_vline(data = aggdf, aes(xintercept = bias), col = "blue", linetype = "dashed") +
    facet_grid(Mixing+Moderation+Separation+Mean_Differences~LatentClass) +
    coord_cartesian(xlim = c(-1,1)) + 
    ggtitle(label = paste0(methods_vec[i],", Small N and No Correlation")) 
  print(p)
}
dev.off()




# Summaries

tmp<-pblapply(X = seq(1,length(parameter_files)),
              FUN = function(x){load(file = paste0(results_wd,"/",summaries_files[x]));return(summaries_df)})
summaries_df = rbindlist(tmp)

KL_min_df = ddply(subset(summaries_df,Imputation_Procedure=="No Missingness"), 
                  .(Sample_Size,Mixing,Moderation,Separation,Mean_Differences,Association), summarize,
              KL_min = weighted.mean(KL, se_KL^(-2)))
summaries_df = left_join(summaries_df, KL_min_df)

summaries_df = transform(summaries_df, d_KL = KL-KL_min)

pdf(file = paste0(stage_wd,"/KL-Density-LargeN-Correlation.pdf"), width = 24, height = 8.8)
methods_vec = unique(means_df$Imputation_Procedure)
  df = subset(summaries_df,Sample_Size == "Large N" & Association == "Correlation")
  aggdf = ddply(df,.(Imputation_Procedure,Mixing,Moderation,Separation,Mean_Differences), summarize, avg_KL = mean(d_KL, na.rm = TRUE))
  p = ggplot(df, aes(d_KL, fill = Imputation_Procedure, col = Imputation_Procedure)) + 
    geom_histogram(bins=50, alpha = 0.5) +
    geom_vline(xintercept = 0, col = "black", size = 1) + 
    geom_vline(data = aggdf, aes(xintercept = avg_KL), col = "red", linetype = "dashed") +
    facet_grid(Imputation_Procedure~Mixing+Moderation+Separation+Mean_Differences) +
    ggtitle(label = paste0("Large N, No Correlation")) + 
    coord_cartesian(xlim = c(0,0.05))
  print(p)
dev.off()


pdf(file = paste0(stage_wd,"/KL-Density-SmallN-Correlation.pdf"), width = 24, height = 8.8)
methods_vec = unique(means_df$Imputation_Procedure)
df = subset(summaries_df,Sample_Size == "Small N" & Association == "Correlation")
aggdf = ddply(df,.(Imputation_Procedure,Mixing,Moderation,Separation,Mean_Differences), summarize, avg_KL = mean(d_KL, na.rm = TRUE))
p = ggplot(df, aes(d_KL, fill = Imputation_Procedure, col = Imputation_Procedure)) + 
  geom_histogram(bins=50, alpha = 0.5) +
  geom_vline(xintercept = 0, col = "black", size = 1) + 
  geom_vline(data = aggdf, aes(xintercept = avg_KL), col = "red", linetype = "dashed") +
  facet_grid(Imputation_Procedure~Mixing+Moderation+Separation+Mean_Differences) +
  ggtitle(label = paste0("Small N, No Correlation")) + 
  coord_cartesian(xlim = c(0,0.2))
print(p)
dev.off()

pdf(file = paste0(stage_wd,"/KL-Density-SmallN-NoCorrelation.pdf"), width = 16, height = 8.8)
methods_vec = unique(means_df$Imputation_Procedure)
df = subset(summaries_df,Sample_Size == "Small N" & Association == "No Correlation")
aggdf = ddply(df,.(Imputation_Procedure,Mixing,Moderation,Separation,Mean_Differences), summarize, avg_KL = mean(d_KL, na.rm = TRUE))
p = ggplot(df, aes(d_KL, fill = Imputation_Procedure, col = Imputation_Procedure)) + 
  geom_histogram(bins=50, alpha = 0.5) +
  geom_vline(xintercept = 0, col = "black", size = 1) + 
  geom_vline(data = aggdf, aes(xintercept = avg_KL), col = "red", linetype = "dashed") +
  facet_grid(Imputation_Procedure~Mixing+Moderation+Separation+Mean_Differences) +
  ggtitle(label = paste0("Small N, No Correlation")) + 
  coord_cartesian(xlim = c(0,0.2))
print(p)
dev.off()
