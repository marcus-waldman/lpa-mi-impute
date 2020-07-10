    rm(list = ls())
    
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
    
    # Directories
    dropbox_wd = "D:/Dropbox"
    #dropbox_wd = "C:/Users/marcu/Dropbox"
    results_wd = paste0(dropbox_wd, "/Dissertation/lpa-mi-impute/stage4c-combine-results")
    stage5b_wd = paste0(dropbox_wd, "/Dissertation/lpa-mi-impute/stage5b-analyze-combined-results")
    environment_wd = paste0(dropbox_wd,"/Dissertation/environmental-variables/")
    
    
    
    # Load in the data
    setwd(results_wd)
    load(file ="parameters-combined-results-lpa-mi-impute.RData")
    load(file ="summaries-combined-results-lpa-mi-impute.RData")
    
    setwd(environment_wd)
    load(file = "labels_for_conditions.Rdata")
    load(file = "data-conditions.Rdata")
    
    # Check missingness rates
    z = 10 # Large sample, no dX, no moderation
    tmp_list<-lapply(X = 1:500, 
                     FUN = function(x){
                       load(paste0("S:/stage0-files/list-observed rep",x," z2 pm1.RData"))
                       return(list_observed$list_obsdf$pm1)
                       })
    big_data = rbindlist(tmp_list)
    rm(tmp_list)
    library(mice)
    md_matrix = md.pattern(big_data[,c("Y1","Y2","Y3","Y4")])
    write.csv(md_matrix, file = paste0(stage5b_wd,"/md-patterns-z2.csv"))
    head(md_matrix)
    R_df = data.frame(md_matrix)
    R_df = R_df[c("Y1","Y2","Y3","Y4")]
    R_df[R_df==1] = Inf
    R_df[R_df==0] = 1
    R_df$Pattern = letters[1:nrow(R_df)]
    R_df = R_df[-nrow(R_df), ]
    names(R_df) = paste0("R_", names(R_df))
    
    R_df$NM = apply(R_df %>% select(starts_with("R_Y")), 1, FUN = function(x){sum(is.finite(x))})
    R_df$NM = mapvalues(R_df$NM, from = 0:3, to = letters[1:4])
    
    big_data = big_data %>% left_join(R_df, by = c("R_Y1","R_Y2","R_Y3","R_Y4"))
    big_data$NM = as.factor(big_data$NM)
    big_data$NM = relevel(big_data$NM, ref = "a")
    
    
    library(nnet)
    test <- multinom(NM ~ Xcom1, data = big_data)
    exp(coef(test)) 
    
    big_data$complete = as.numeric(big_data$NM != "a")
    
    obj_glm = glm(complete~Xcom1, data = big_data, family = "binomial")
    exp(coef(obj_glm))
    setwd(stage5b_wd)
    
    # Merge the environmental conditions
    parameters_combined_df = parameters_combined_df %>% left_join(conditions_labels,by="z")
    summaries_combined_df = summaries_combined_df %>% left_join(conditions_labels,by="z")
    
    # Descriptive statistics
    with(summaries_combined_df,table(z,Imputation_Procedure)) 
      #See that if no FIML then no Pattern Mixture 
      #Makes sense because didn't fit models with no FIML
    
    # Check that I sampled all of the correct conditions. Short answer: No.
    tmp_df = conditions_labels %>% left_join(data_conditions)
    with(tmp_df, table(Mixing,class_size))
    with(tmp_df, table(Sample.Size, N)) #Unequal!
    with(tmp_df, table(Separation,MD))
    with(tmp_df, table(AV.Correlation, rho_YX)) # Many more
    with(tmp_df, table(AV.Mean.Diff, dX))
    with(tmp_df, table(AV.Moderation, C_modifies_YX))
    tmp_df = transform(tmp_df, txt = as.character(paste(Mixing, Separation, AV.Correlation, AV.Mean.Diff, AV.Moderation, sep = ", ")))
    txt_small = subset(tmp_df, N == 300)$txt
    txt_large = subset(tmp_df, N == 1200)$txt
    setdiff(txt_small, txt_large)
        # [1] "Equal Mix., Weakly Separation, AV-No Corr., AV-Mean Class Diff., AV-No Moderation"   
        # [2] "Equal Mix., Strongly Separated, AV-No Corr., AV-Mean Class Diff., AV-No Moderation"  
        # [3] "Equal Mix., Strongly Separated, AV-No Corr., AV-Mean Class Diff., AV-Moderation"     
        # [4] "Equal Mix., Weakly Separation, AV-No Corr., AV-Mean Class Diff., AV-Moderation"      
        # [5] "Unequal Mix., Strongly Separated, AV-No Corr., AV-Mean Class Diff., AV-No Moderation"
        # [6] "Unequal Mix., Weakly Separation, AV-No Corr., AV-Mean Class Diff., AV-No Moderation" 
        # [7] "Unequal Mix., Weakly Separation, AV-No Corr., AV-Mean Class Diff., AV-Moderation"    
        # [8] "Unequal Mix., Strongly Separated, AV-No Corr., AV-Mean Class Diff., AV-Moderation" 
    
    #Suggests that I did the following improper things:
        # 1) In the large sample, I inappropriately did not consider the situation where there were there was no correlation but there did exist class mean differences 
        #     this would result in 2(Mixing)*2(Separation)*1(AV-No Moderation) = 4 additional conditions that I should have considered.
        # 2) In the small sample, I considered redundant conditions in that Both AV-No Corr is a condition which makes it not matter where there is a moderation or not. 
        #   However, both moderation prcoesses are present (e.g., AV-No Moderation and AV-Moderation). 
        #   This means that I considered 8/2=4 additional conditions in the small sample size that I should not have. 
    #*Upshot: As a result ofthis, I am am going to get rid of any conditions with No correlation and not analyze it. 
    
    # Summaries 
    parameters_combined_df = subset(parameters_combined_df, subset = AV.Correlation=="AV-Strong Corr.")
    summaries_combined_df = subset(summaries_combined_df, subset = AV.Correlation=="AV-Strong Corr.")
    
    # Generate value labels for the imputation procedure
    imputation_procedure_table = data.frame(Imputation_Procedure = as.character(c("No Missingness","FIML","Pattern Mixture","Amelia","PMM","CART","RF")))
    
    imputation_procedure_table = transform(imputation_procedure_table, 
                                           pno = 0:6)
    imputation_procedure_table = transform(imputation_procedure_table,
                                           Strategy = ordered(pno, levels = c(0:6), label = c("Complete data","FIML","Pattern Mixture Model","EMS-MVN.","MICE-PMM","MICE-CART","MICE-RF")))
    imputation_procedure_table$Imputation_Procedure = as.character(imputation_procedure_table$Imputation_Procedure)                                       
    
    
    parameters_combined_df = parameters_combined_df %>% left_join(imputation_procedure_table, by = "Imputation_Procedure")
    summaries_combined_df = summaries_combined_df %>% left_join(imputation_procedure_table, by = "Imputation_Procedure")
    
    with(parameters_combined_df, table(Strategy, Imputation_Procedure))
    with(parameters_combined_df, table(Strategy, Imputation_Procedure))
    
#pdf(file = paste0(stage5b_wd,"/binder-of-plots-17-june-2019.pdf"), width = 7.5, height = 10)    
    
    # Marginal results
    marginal_KL = ddply(summaries_combined_df, .(Strategy), summarise, 
                        KL = weighted.mean(KL, w = se_KL^(-2), na.rm = TRUE))
    print(marginal_KL)
    
    
    
    ##########################################################################
    # One-way results
    ##########################################################################
    oneway_KL_samplesize = ddply(summaries_combined_df, .(Sample.Size,Strategy), summarise, 
                                 KL_est = weighted.mean(KL, w = se_KL^(-2), na.rm = TRUE), 
                                 KL_ub = quantile(KL, probs = 0.750, na.rm = TRUE), 
                                 KL_lb = quantile(KL, probs = 0.250, na.rm = TRUE))
    print(oneway_KL_samplesize)
    
    tmp1_df = subset(oneway_KL_samplesize, Strategy == "Complete data") %>% select("Sample.Size","KL_est")
    names(tmp1_df)[2] = c("KL_complete")
    
    tmp2_df = subset(oneway_KL_samplesize, Strategy == "FIML") %>% select("Sample.Size","KL_est")
    names(tmp2_df)[2] = c("KL_FIML")
    
    oneway_KL_samplesize = oneway_KL_samplesize %>% left_join(tmp1_df) %>% left_join(tmp2_df) %>%
      transform(Percent.Change.KL = 100*(KL_est-KL_FIML)/(KL_FIML-KL_complete)) %>% 
      transform(Percent.Change.KL.UB =  100*(KL_ub-KL_FIML)/(KL_FIML-KL_complete)) %>%
      transform(Percent.Change.KL.LB =  100*(KL_lb-KL_FIML)/(KL_FIML-KL_complete))
      
    # 1A) Simulation 1 Imputation:
    plot_1a = ggplot(subset(oneway_KL_samplesize, Strategy != "Complete data" & Strategy != "Pattern Mixture Model"), 
                     aes(x = Strategy, y = -1*Percent.Change.KL, ymax = -1*Percent.Change.KL.UB, ymin = -1*Percent.Change.KL.LB))
    plot_1a = plot_1a + geom_abline(slope = 0, intercept = 100, col = "gray50", linetype = 1, size = 3) + geom_crossbar(fatten = 3, col = "black") 
    plot_1a = plot_1a + labs(x = "", y = "",title = "% Reduction in Kullback-Leibler Divergence")
    plot_1a = plot_1a + coord_cartesian(ylim = c(-75,100))
    plot_1a = plot_1a + facet_grid(~Sample.Size)
    plot_1a = plot_1a + theme_bw()
    plot_1a = plot_1a + theme(legend.position = "top", plot.title = element_text(hjust = 0.5), 
                                legend.title = element_blank(), axis.text.x = element_text(angle=90, hjust=1))
    print(plot_1a) 
    ggsave(filename = "oneway-pct-reduc-kl.png", plot = plot_1a, width = 5, height = 5)
    
    # # Two-way results
    # twoway_KL_samplesize_mixing = ddply(summaries_combined_df, .(Sample.Size, Mixing, Strategy), summarise, 
    #                                     KL_est = weighted.mean(KL, w = se_KL^(-2), na.rm = TRUE), 
    #                                     KL_ub = quantile(KL, probs = 0.750, na.rm = TRUE), 
    #                                     KL_lb = quantile(KL, probs = 0.250, na.rm = TRUE))
    # 
    # tmp1_df = subset(twoway_KL_samplesize_mixing, Strategy == "Complete data") %>% select("Sample.Size","Mixing","KL_est")
    # names(tmp1_df)[3] = c("KL_complete")
    # 
    # tmp2_df = subset(twoway_KL_samplesize_mixing, Strategy == "FIML") %>% select("Sample.Size","Mixing","KL_est")
    # names(tmp2_df)[3] = c("KL_FIML")          
    #           
    # twoway_KL_samplesize_mixing = twoway_KL_samplesize_mixing %>% left_join(tmp1_df) %>% left_join(tmp2_df) %>%
    #   transform(Percent.Change.KL = 100*(KL_est-KL_FIML)/(KL_FIML-KL_complete)) %>% 
    #   transform(Percent.Change.KL.UB =  100*(KL_ub-KL_FIML)/(KL_FIML-KL_complete)) %>%
    #   transform(Percent.Change.KL.LB =  100*(KL_lb-KL_FIML)/(KL_FIML-KL_complete))
    # 
    # # 1B) Simulation 1 Imputation:
    # plot_1b = ggplot(subset(twoway_KL_samplesize_mixing, Strategy != "Complete data" & Strategy != "Pattern Mixture Model"), 
    #                  aes(x = Strategy, y = -1*Percent.Change.KL, ymax = -1*Percent.Change.KL.UB, ymin = -1*Percent.Change.KL.LB, 
    #                      fill = Mixing, linetype = Mixing))
    # plot_1b = plot_1b + geom_abline(slope = 0, intercept = 100, col = "gray50", linetype = 1, size = 3) 
    # plot_1b = plot_1b + geom_crossbar(fatten = 2, alpha = 0.25)
    # plot_1b = plot_1b + labs(x = "", y = "",title = "% Reduction in Kullback-Leibler Divergence")
    # plot_1b = plot_1b + coord_cartesian(ylim = c(-75,100))
    # plot_1b = plot_1b + facet_grid(~Sample.Size)
    # plot_1b = plot_1b + theme_bw()
    # print(plot_1b) 
    
    
    
    
    # Three-way results
    threeway_KL_samplesize_mixing_sep = ddply(summaries_combined_df, .(Sample.Size, Mixing, Separation, Strategy), summarise, 
                                        KL_est = weighted.mean(KL, w = se_KL^(-2), na.rm = TRUE), 
                                        KL_ub = quantile(KL, probs = 0.750, na.rm = TRUE), 
                                        KL_lb = quantile(KL, probs = 0.250, na.rm = TRUE), 
                                        ARIU = mean(ARIU, na.rm = TRUE))
    
    
    
    tmp1_df = subset(threeway_KL_samplesize_mixing_sep, Strategy == "Complete data") %>% select("Sample.Size","Mixing","Separation","KL_est")
    names(tmp1_df)[4] = c("KL_complete")
    
    tmp2_df = subset(threeway_KL_samplesize_mixing_sep, Strategy == "FIML") %>% select("Sample.Size","Mixing","Separation","KL_est")
    names(tmp2_df)[4] = c("KL_FIML")          
    
    threeway_KL_samplesize_mixing_sep = threeway_KL_samplesize_mixing_sep %>% left_join(tmp1_df) %>% left_join(tmp2_df) %>%
      transform(Percent.Change.KL = 100*(KL_est-KL_FIML)/(KL_FIML-KL_complete)) %>% 
      transform(Percent.Change.KL.UB =  100*(KL_ub-KL_FIML)/(KL_FIML-KL_complete)) %>%
      transform(Percent.Change.KL.LB =  100*(KL_lb-KL_FIML)/(KL_FIML-KL_complete))
    
    
    
    # 1C) Simulation 1 Imputation:
    plot_1c = ggplot(subset(threeway_KL_samplesize_mixing_sep, Strategy != "Complete data" & Strategy != "Pattern Mixture Model"), 
                     aes(x = Strategy, y = -1*Percent.Change.KL, ymax = -1*Percent.Change.KL.UB, ymin = -1*Percent.Change.KL.LB, 
                         fill = Mixing, linetype = Mixing))
    plot_1c = plot_1c + geom_abline(slope = 0, intercept = 100, col = "gray50", linetype = 1, size = 3) 
    plot_1c = plot_1c + geom_crossbar(fatten = 3, alpha = 0.1)
    plot_1c = plot_1c + labs(x = "", y = "",title = "% Reduction in Kullback-Leibler Divergence")
    plot_1c = plot_1c + coord_cartesian(ylim = c(-75,100))
    plot_1c = plot_1c + facet_grid(Sample.Size~Separation)
    plot_1c = plot_1c + theme_bw() 
    plot_1c = plot_1c + theme(legend.position = "top", plot.title = element_text(hjust = 0.5), 
                              legend.title = element_blank(), axis.text.x = element_text(angle=90, hjust=1))
    print(plot_1c) 
    ggsave(filename = "Pct-Reduc-KL.png", plot = plot_1c, device = "png", width = 5, height = 5)
    
    
    
    
table_at1_kl = threeway_KL_samplesize_mixing_sep %>% 
  filter(Strategy != "Complete data") %>% 
  select("Sample.Size","Mixing","Separation","Strategy","Percent.Change.KL") %>% 
  spread(Sample.Size,Percent.Change.KL) 
table_at1_kl = table_at1_kl %>% arrange(desc(Mixing), Separation,Strategy)
names(table_at1_kl)[seq(ncol(table_at1_kl)-1, ncol(table_at1_kl))] = c("KL.Small.N","KL.Large.N")
table_at1_kl$KL.Large.N = -1*table_at1_kl$KL.Large.N
table_at1_kl$KL.Small.N = -1*table_at1_kl$KL.Small.N
write.csv(table_at1_kl, file = "kl-results-for-table-at1.csv", row.names = FALSE)
    
    
    plot_1d = ggplot(subset(threeway_KL_samplesize_mixing_sep, Strategy != "Complete data" & Strategy != "Pattern Mixture Model" & Strategy != "FIML"), 
                     aes(x = Strategy, y = ARIU, shape = Mixing))
    plot_1d = plot_1d + geom_point(size = 3) 
    plot_1d = plot_1d + labs(x = "", y = "",title = "Average Relative Increase in Variance")
    plot_1d = plot_1d + facet_wrap(~Sample.Size+Separation, scales = "free_y")
    plot_1d = plot_1d + theme_bw() 
    plot_1d = plot_1d + theme(legend.position = "top", plot.title = element_text(hjust = 0.5), 
                              legend.title = element_blank(), axis.text.x = element_text(angle=90, hjust=1))
    print(plot_1d) 
    

    # Bias in means
    means_Y_df = ddply(parameters_combined_df %>% filter(paramHeader=="Means"&startsWith(param,"Y")), 
                       .(Sample.Size, Mixing,Separation,LatentClass,param,Strategy), summarise, 
                    bias = weighted.mean(deviance, na.rm = TRUE))
    
    means_Y_FIML = means_Y_df %>% filter(Strategy == "FIML") %>% select("Sample.Size", "Mixing", "Separation", "LatentClass", "param","bias")
    names(means_Y_FIML)[names(means_Y_FIML) == "bias"] = "bias_FIML" 
    
    means_Y_COMPLETE = means_Y_df %>% filter(Strategy == "Complete data") %>% select("Sample.Size", "Mixing", "Separation", "LatentClass", "param","bias")
    names(means_Y_COMPLETE)[names(means_Y_COMPLETE) == "bias"] = "bias_COMPLETE" 
    
    means_Y_df = means_Y_df %>% left_join(means_Y_FIML) %>%
      left_join(means_Y_COMPLETE) %>% transform(rr_bias = bias)
    
    for ( tmp.Sample.Size in c("Small Sample (N=300)", "Large Sample (N=1200)") ){
      for (tmp.Mixing in c("Equal Mix.", "Unequal Mix.")){
        
        df_means_tmp = means_Y_df %>% filter(Strategy != "FIML", Strategy != "Pattern Mixture Model", 
                                             Sample.Size == tmp.Sample.Size, Mixing == tmp.Mixing)
        p_tmp = ggplot(df_means_tmp, 
               aes(x = LatentClass, y = abs(rr_bias), col = Separation, fill = Separation)) + 
          geom_hline(yintercept = 0, size = 1, alpha =1) + 
          #geom_point(col = "black", shape = 0) + 
          facet_grid (param~Strategy) + 
          geom_segment(data = means_Y_df %>% filter(Strategy != "FIML", Strategy != "Pattern Mixture Model", Sample.Size == tmp.Sample.Size, Mixing == tmp.Mixing),
                    aes(x = LatentClass, y = abs(bias_FIML), xend = LatentClass, yend = abs(rr_bias)), 
                    arrow = arrow(length = unit(0.075, "npc"), type = "closed"), size = 0.75)  +  
          geom_point(data = means_Y_df %>% filter(Strategy != "FIML", Strategy != "Pattern Mixture Model",  Sample.Size == tmp.Sample.Size, Mixing == tmp.Mixing), 
                     aes(x = LatentClass, y = abs(bias_FIML)), shape = 21, size = 2) + 
          coord_cartesian(ylim = c(-0.02,0.15)) + 
          labs(title = paste0(tmp.Sample.Size, " & ", tmp.Mixing), x = "Class", y = "|Bias|") + 
          scale_fill_hue(l=40, c=35) + 
          scale_color_hue(l=40, c=35) +
          theme(legend.position = "top", plot.title = element_text(hjust = 0.5))
        
        ggsave(filename = paste0("AbsBias-",tmp.Sample.Size,"-",tmp.Mixing,"png"), plot = p_tmp, device = "png", width = 10, height = 7.5)
      
        
      }
    }
    
    means_Y_df = transform(means_Y_df, rr_bias = bias/bias_FIML)
    table_Y_bias = means_Y_df %>% select(Sample.Size, Mixing, Separation, Strategy,LatentClass, param,rr_bias) %>% 
      transform(LatentClass = paste("Class", LatentClass)) %>% 
      group_by(Sample.Size,Mixing,Separation,Strategy) %>% 
      spread(key = LatentClass, value = rr_bias) %>%
      arrange(param, Sample.Size, Mixing, Separation, Strategy)
    write.csv(table_Y_bias, "table_Y_bias.csv", row.names = FALSE)
    
    
    ids_hash1 = which(parameters_combined_df$param == "C#1")
    ids_hash2 = which(parameters_combined_df$param == "C#2")
    hash_df = data.frame(id = 1:length(ids_hash1),
                         idx_1 = ids_hash1, 
                         idx_2 = ids_hash2) %>% 
             data.table() %>% 
             melt.data.table(id.vars = "id", measure.vars = c("idx_1","idx_2"), value.name = "ids_hash") %>% 
             arrange(id, variable) %>% 
             transform(gamma = parameters_combined_df$est[ids_hash], 
                       value = parameters_combined_df$value[ids_hash], 
                       rep = parameters_combined_df$rep[ids_hash], 
                       z = parameters_combined_df$z[ids_hash], 
                       pm = parameters_combined_df$pm[ids_hash], 
                       pva = parameters_combined_df$pva[ids_hash], 
                       kfit = parameters_combined_df$kfit[ids_hash], 
                       MM = parameters_combined_df$MM[ids_hash], 
                       data_type = parameters_combined_df$data_type[ids_hash], 
                       Imputation_Procedure = parameters_combined_df$Imputation_Procedure[ids_hash], 
                       Mixing = parameters_combined_df$Mixing[ids_hash], 
                       Sample.Size = parameters_combined_df$Sample.Size[ids_hash], 
                       Separation = parameters_combined_df$Separation[ids_hash], 
                       AV.Correlation = parameters_combined_df$AV.Correlation[ids_hash], 
                       AV.Mean.Diff = parameters_combined_df$AV.Mean.Diff[ids_hash],
                       AV.Moderation = parameters_combined_df$AV.Moderation[ids_hash], 
                       pno = parameters_combined_df$pno[ids_hash], 
                       Strategy = parameters_combined_df$Strategy[ids_hash])
    
    # library(parallel)
    # detectCores()
    # cl = makeCluster(detectCores())
    # clusterExport(cl, varlist = c("hash_df"))
    # tmp_list <- pblapply(X = unique(hash_df$id), 
    #                    FUN = function(x){
    #                           library(lpa.mi.src); library(tidyverse); library(data.table)
    #                           hash_x = hash_df %>% filter(id==x)
    #                           pi_x = with(hash_x, lpa.mi.src::gamma2pi(gamma))
    #                           value_x = with(hash_x, lpa.mi.src::gamma2pi(value))
    #                           pi_df = data.table(value = value_x, est = pi_x, id = hash_x$id[1], LatentClass = 1:length(value_x)) %>% 
    #                             transform(deviance = est-value) %>% 
    #                             left_join(hash_x[1,c("id","rep","z","pm","pva","kfit","MM","data_type","Imputation_Procedure","Mixing","Sample.Size", 
    #                                                  "Separation","AV.Correlation","AV.Mean.Diff","AV.Moderation","pno","Strategy")], 
    #                                       by = "id")
    #                           return(pi_df)
    #                 }, 
    #                 cl = cl)
    # stopCluster(cl)
    # pi_df = rbindlist(tmp_list)
    # rm(tmp_list)
    #save(pi_df, file = "pi_df.RData")
    load("pi_df.RData")
    
    
    means_PI_df = ddply(pi_df, 
                       .(Sample.Size, Mixing,Separation,LatentClass,Strategy), summarise, 
                       bias = mean(deviance, na.rm = TRUE))
    
    means_PI_FIML = means_PI_df %>% filter(Strategy == "FIML") %>% select("Sample.Size", "Mixing", "Separation","LatentClass","bias")
    names(means_PI_FIML)[names(means_PI_FIML) == "bias"] = "bias_FIML" 
    
    means_PI_COMPLETE = means_PI_df %>% filter(Strategy == "Complete data") %>% select("Sample.Size", "Mixing", "Separation", "LatentClass","bias")
    names(means_PI_COMPLETE)[names(means_PI_COMPLETE) == "bias"] = "bias_COMPLETE" 
    
    means_PI_df = means_PI_df %>% left_join(means_PI_FIML) %>%
      left_join(means_PI_COMPLETE) %>% transform(rr_bias = bias)
    
    
    for ( tmp.Sample.Size in c("Small Sample (N=300)", "Large Sample (N=1200)") ){
      for (tmp.Mixing in c("Equal Mix.", "Unequal Mix.")){
#tmp.Sample.Size = "Large Sample (N=1200)"
#tmp.Mixing = "Equal Mix."
        df_means_tmp = means_PI_df %>% filter(Strategy != "FIML", Sample.Size == tmp.Sample.Size, Mixing == tmp.Mixing)
        ggplot(df_means_tmp, 
               aes(x = as.ordered(LatentClass), y = abs(rr_bias), col = Separation, fill = Separation)) + 
          geom_hline(yintercept = 0, size = 1, alpha =1) + 
          #geom_point(col = "black", shape = 0) + 
          facet_grid (Separation~Strategy) + 
          geom_segment(data = means_PI_df %>% filter(Strategy != "FIML", Sample.Size == tmp.Sample.Size, Mixing == tmp.Mixing),
                       aes(x = as.ordered(LatentClass), y = abs(bias_FIML), xend = as.ordered(LatentClass), yend = abs(rr_bias)), 
                       arrow = arrow(length = unit(0.075, "npc"), type = "closed"), size = 0.75)  +  
          geom_point(data = means_PI_df %>% filter(Strategy != "FIML", Sample.Size == tmp.Sample.Size, Mixing == tmp.Mixing), 
                     aes(x = as.ordered(LatentClass), y = abs(bias_FIML)), shape = 21, size = 2) + 
          coord_cartesian(ylim = c(0,0.1)) + 
          labs(title = paste0(tmp.Sample.Size, " & ", tmp.Mixing), x = "Class", y = "|Bias| in Mix. Proportions") + 
          scale_fill_hue(l=40, c=35) + 
          scale_color_hue(l=40, c=35) +
          theme(legend.position = "bottom", plot.title = element_text(hjust = 0.5))
        ggsave(paste0("pi-AbsBias-",tmp.Sample.Size,"-",tmp.Mixing,"png"), device = "png", width = 10, height = 5.5)
      }
    }
    
    means_PI_df = transform(means_PI_df, rr_bias = bias/bias_FIML)
    table_PI_bias = means_PI_df %>% select(Sample.Size, Mixing, Separation, Strategy,LatentClass, rr_bias) %>% 
      transform(LatentClass = paste("Class", LatentClass)) %>% 
      group_by(Sample.Size,Mixing,Separation,Strategy) %>% 
      spread(key = LatentClass, value = rr_bias)
    write.csv(table_PI_bias, "table_PI_bias.csv", row.names = FALSE)
    
    
    
    
    hi = parameters_combined_df %>% 
      filter(paramHeader=="Means" & startsWith(param,"Y") & Strategy!="Pattern Mixture Model") %>%
      select(Sample.Size, Mixing,Separation,param,Strategy,deviance,covered) %>%
      group_by(Sample.Size, Mixing,Separation,Strategy,param) %>% 
      summarize_all(funs(mean))
    
    df_1e = parameters_combined_df %>% 
              filter(paramHeader=="Means" & startsWith(param,"Y") & Strategy!="Pattern Mixture Model") %>%
              select(Sample.Size,Strategy,LatentClass,param,deviance,covered) %>%
              group_by(Sample.Size,Strategy,LatentClass,param) %>% 
              summarize_all(funs(mean)) %>%
              transform(Accept.Cov = as.factor(ifelse(covered>=0.925&covered<=0.975, "Yes","No")))
    
    plot_1e = ggplot(df_1e, aes(x = as.ordered(param), y = covered, label = LatentClass, col = Accept.Cov))
    #plot_1e = plot_1e + geom_abline(slope = 0, intercept  = 0.95)
    plot_1e = plot_1e + geom_abline(slope = 0, intercept  = 0.975,linetype=3)
    plot_1e = plot_1e + geom_abline(slope = 0, intercept  = 0.925,linetype=3)
    plot_1e = plot_1e + geom_text(size = 4)
    plot_1e = plot_1e + facet_grid(Strategy~Sample.Size, scales = "free_y") 
    plot_1e = plot_1e + guides(color = "none")
    plot_1e = plot_1e + labs(title = "95% CI Coverage Rate", x = "", y = "")
    plot_1e = plot_1e + theme_bw()
    plot_1e = plot_1e + theme(plot.title = element_text(hjust = 0.5))
    plot_1e = plot_1e + geom_hline(yintercept = 1)
    print(plot_1e)
    ggsave(filename = "CI95.png", plot = plot_1e, width = 5, height = 7)
    
    
    df_1f = parameters_combined_df %>% 
      filter(paramHeader=="Means" & startsWith(param,"Y") & Strategy!="Pattern Mixture Model") %>%
      select(Sample.Size,Mixing,Separation,Strategy,LatentClass,param,deviance,covered) %>%
      group_by(Sample.Size,Mixing,Separation,Strategy,LatentClass,param) %>% 
      summarize_all(funs(mean)) %>%
      mutate(Accept.Cov = ifelse(covered>=0.925&covered<=0.975, "Yes","No"))
    plot_1f = ggplot(df_1f, aes(x = as.ordered(param), y = covered, 
                                label = LatentClass, fill = Mixing, color = Mixing))
    plot_1f = plot_1f + geom_abline(slope = 0, intercept  = 0.95)
    plot_1f = plot_1f + geom_text()
    plot_1f = plot_1f + facet_grid(Sample.Size+Separation~Strategy+Mixing) 
    print(plot_1f)
    
    
    
    
    # Bias in mixing proportions
    
    
    df_1g = parameters_combined_df %>% 
      filter(paramHeader=="Means" & startsWith(param,"C#") & Strategy!="Pattern Mixture Model") %>%
      select(Sample.Size,Strategy,Mixing,Separation,param,deviance,covered) %>%
      group_by(Sample.Size,Strategy,Separation,Mixing,param) %>% 
      summarize_all(funs(mean)) %>%
      mutate(Accept.Cov = ifelse(covered>=0.925&covered<=0.975, "Yes","No"))
    
    plot_1g= ggplot(df_1g, aes(x = Strategy, y = covered, label = param, hjust = param, col = Accept.Cov, shape = Mixing))
    plot_1g = plot_1g + geom_abline(slope = 0, intercept  = 0.975,linetype=3)
    plot_1g = plot_1g + geom_abline(slope = 0, intercept  = 0.925,linetype=3)
    plot_1g = plot_1g + geom_point(size = 3)
    plot_1g = plot_1g + facet_grid(Sample.Size~Separation+param) 
    plot_1g = plot_1g + labs(title = "Log-odds Ratios: 95% CI Coverage Rate", x = "", y = "")
    plot_1g = plot_1g + theme_bw()
    plot_1g = plot_1g + theme(legend.position = "top", plot.title = element_text(hjust = 0.5), 
                              legend.title = element_blank(), axis.text.x = element_text(angle=90, hjust=1))
    print(plot_1g)
    
# dev.off()   

    
    # parameters_combined_df %>% 
    #   filter(paramHeader=="Means" & startsWith(param,"C#") & Strategy!="Pattern Mixture Model") %>%
    #   select(Sample.Size, Mixing,Separation,param,Strategy,deviance,covered) %>%
    #   group_by(Sample.Size, Mixing,Separation,Strategy,param) %>% 
    #   summarize_all(funs(mean))
    # 

df_1f_v2 = parameters_combined_df %>% 
  filter(paramHeader=="Means" & startsWith(param,"Y")) %>%
  select(Sample.Size,Mixing,Separation,Strategy,LatentClass,param,deviance,covered) %>%
  group_by(Sample.Size,Mixing,Separation,Strategy,LatentClass,param) %>% 
  summarize_all(funs(mean)) %>%
  mutate(Accept.Cov = ifelse(covered>=0.925&covered<=0.975, "Yes","No"))
    
table_at1_cov_means = ddply(data.frame(df_1f_v2), .(Sample.Size,Mixing,Separation,Strategy), summarise, 
                      Avg.Cov.Mean = mean(covered, na.rm = TRUE)) %>% 
                      spread(Sample.Size, Avg.Cov.Mean)
names(table_at1_cov_means)[seq(ncol(table_at1_cov_means)-1,ncol(table_at1_cov_means))] = c("Small.Sample.Cov.Means", "Large.Sample.Cov.Means")
                      

df_1g_v2 = parameters_combined_df %>% 
  filter(paramHeader=="Means" & startsWith(param,"C#")) %>%
  select(Sample.Size,Strategy,Mixing,Separation,param,deviance,covered) %>%
  group_by(Sample.Size,Strategy,Separation,Mixing,param) %>% 
  summarize_all(funs(mean)) %>%
  mutate(Accept.Cov = ifelse(covered>=0.925&covered<=0.975, "Yes","No"))



table_at1_cov_lodds = ddply(data.frame(df_1g_v2), .(Sample.Size,Mixing,Separation,Strategy), summarise, 
                            Avg.Cov.Lodd = mean(covered, na.rm = TRUE)) %>% 
                      spread(Sample.Size, Avg.Cov.Lodd)
names(table_at1_cov_lodds)[seq(ncol(table_at1_cov_lodds)-1,ncol(table_at1_cov_lodds))] = c("Small.Sample.Cov.Lodds", "Large.Sample.Cov.Lodds")

table_at1_cov = table_at1_cov_means %>% left_join(table_at1_cov_lodds)

#write.csv(table_at1_cov, file = "coverage-results-for-table-at1.csv", row.names = FALSE)



