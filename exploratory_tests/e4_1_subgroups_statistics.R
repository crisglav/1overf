# Compare age-corrected aperiodic exponents between healthy participants and the two
# biggests subgroups of patients (Chronic back pain patients = CBP) and
# chronic widespread pain patients (CWP)
#
# Cristina Gil, Flaminia Palloti, TUM, 27.02.2024

rm(list=ls())
library(BayesFactor)
library(dplyr)
library(effsize)

# Load data from the main analysis' results.
file <- "/rechenmagd3/Experiments/2023_1overf/results/statistics/exp_PFC_real.csv"
participants = read.csv(file, header = TRUE, sep = ",")

# Output data
file_cwp <- "/rechenmagd3/Experiments/2023_1overf/results/statistics/e4_cwp.csv"
file_cbp <- "/rechenmagd3/Experiments/2023_1overf/results/statistics/e4_cbp.csv"

##################################
# TRANSFORMATIONS TO THE DATA
##################################
# Replace all Nan and n/a values with NA 
participants <- participants %>%
  mutate_all(~ifelse(. %in% c("NaN","n/a"), NA, .))

# Create a new variable avg_pain_new that takes the current pain for all the patients that don't have a value in avg_pain
participants <- participants %>%
  mutate(avg_pain_new = ifelse(is.na(avg_pain),curr_pain, avg_pain))

# Code gender as a binary variable
participants <- participants %>%
  mutate(isfemale = ifelse(gender == "f", 1,0))

# # Regress out age from aperiodic exponents and code it as a new variable
# model <- lm(participants$exp_PFC ~ participants$age)
# participants$res_apexp <- resid(model)
# 
# # Regress out age from pain ratings and code it as a new variable
# model <- lm(participants$avg_pain_new ~ participants$age)
# participants <- participants %>%
#   mutate(res_pain = ifelse(group == "pa",resid(model),NA))

#################################

# Separate participants by group and diagnosis
hc = subset(participants, group == 'hc')
cbp = subset(participants, diagnosis == 'CBP')
cwp = subset(participants, diagnosis == 'CWP')

# Auxiliary function that performs a bayesian ttest between groups in 
# age and bender and returns the correspondent BF
compute_bf <- function (patients, healthy){
  bttest_age <- ttestBF(patients$age, healthy$age, mu = 0, paired = FALSE, rscale = "medium", posterior = FALSE)
  bf_age <- as.vector(bttest_age)
  bttest_gender <- ttestBF(patients$isfemale, healthy$isfemale, mu = 0, paired = FALSE, rscale = "medium", posterior = FALSE)
  bf_gender <- as.vector(bttest_gender)
  return (list(bf_age = bf_age, bf_gender = bf_gender))
}


############################################################
# CWP vs Healthy subgroup 
############################################################
if (!file.exists(file_cwp)){
  # CWP patients from FM study (Tiemann et. al, 2012), original ID is "sub-FMpaXX"
  cwp_tiemann <- cwp[grepl("sub-FMpa[0-9]{2}",cwp$original_id),]
  hc_tiemann <- hc[grepl("sub-FMhc[0-9]{2}",hc$original_id),]
  
  # CWP patients from NCCP study (Ta Dinh et. al, 2019), original ID is "sub-NCCPpaXX"
  cwp_tadinh <- cwp[grepl("sub-NCCPpa[0-9]{2}",cwp$original_id),]
  hc_tadinh <- hc[grepl("sub-NCCPhc[0-9]{2}",hc$original_id),]
  
  # CWP patients from longitudinal study (Heitmann et. al, 2022), original ID is "sub-XX".
  # This dataset does not have a matching healthy group but it was recorded with the same settings as Ta Dinh dataset
  cwp_heitmann <- cwp[grepl("sub-[0-9]{2}",cwp$original_id),]
  
  # Select a subset of healthy participants that match the characteristics of the CWP participants
  # i.e. that belonged to the same initial project and that don't differ in age or gender
  # We follow an iterative procedure that selects a random subset of healthy participants 
  # until we get evidence that there's no difference between groups in age and gender (BF<0.33)
  
  bf_age <- Inf
  bf_gender <- Inf
  selected_controls <- NULL
  counter <- 0
  while(bf_age >= 0.33 | bf_gender >= 0.33 & counter < 100){
    
    # select subsample from healthy population of size of the corresp. patient population
    hc_cwp <- rbind(
      hc_tiemann[sample(nrow(hc_tiemann), nrow(cwp_tiemann), replace = FALSE),],  
      hc_tadinh[sample(nrow(hc_tadinh), nrow(cwp_tadinh)+nrow(cwp_heitmann), replace = FALSE),]) 
    
    bf_results <- compute_bf(cwp, hc_cwp)
    
    bf_age <- bf_results$bf_age
    bf_gender <- bf_results$bf_gender
    counter <- counter +1
  }
  print("Exit Loop")
  
  data_cwp <- rbind(hc_cwp,cwp)
  
  # Regress out age from the aperiodic exponents in the subgroup of participants
  model <- lm(data_cwp$exp_PFC ~ data_cwp$age)
  res_apexp <- resid(model)
  res_cwp <-res_apexp[data_cwp$group == "pa"]
  res_hc_cwp <-res_apexp[data_cwp$group == "hc"]
  data_cwp$res_apexp_h1 = res_apexp
  
  # Regress out age from the aperiodic exponents in the cwp patients
  cwp <- data_cwp[data_cwp$group == "pa",]
  model <- lm(cwp$exp_PFC ~ cwp$age)
  res_apexp_cwp <-resid(model)
  # Regress out age from pain ratings in the cwp patients
  model <- lm(cwp$avg_pain_new ~ cwp$age)
  res_pain_cwp <-resid(model)
  data_cwp$res_apexp_h2[data_cwp$group == "pa"] = res_apexp_cwp
  data_cwp$res_pain_h2[data_cwp$group == "pa"] = res_pain_cwp
  
  write.table(data_cwp, file = file_cwp, sep=",", row.names = FALSE)
  
  
} else{
  data_cwp <- read.csv(file_cwp, header = TRUE, sep = ",")
  # hc_cwp <- data_cwp[data_cwp$group == "hc",]
  # cwp <- data_cwp[data_cwp$group == "pa",]
}


# H1. Bayesian two sided independent samples t-test between groups on the residuals
bttest = ttestBF(res_hc_cwp, res_cwp, mu = 0, paired = FALSE, rscale = "medium", posterior = FALSE)
post = posterior(bttest,iterations = 1000)
bf_cwp_h1 = as.data.frame(bttest)$bf
postdelta_cwp_h1 = median(post[,"delta"])
# Frequentist t-test and effect size (Cohen's d)
fttest = t.test(res_hc_cwp, res_cwp, alternative = "two.sided", var.equal = FALSE)
pvalue_cwp_h1 = fttest$p.value
eff = cohen.d(res_hc_cwp,res_cwp)
cohens_d_cwp_h1 = eff$estimate

# H2. Bayesian correlation between age-corrected aperiodic exponents and age-corrected pain ratings
# in the CWP dataset
bcorrelation = correlationBF(y = res_apexp_cwp, x = res_pain_cwp, rscale = "ultrawide")
post = posterior(bcorrelation,iterations = 1000)
bf_cwp_h2 = as.data.frame(bcorrelation)$bf
postrho_cwp_h2 = median(post[,"rho"])
# Pearson's correlation
corr = cor.test(y = res_apexp_cwp, x = res_pain_cwp)
r_cwp_h2 = as.data.frame(corr$estimate)$cor
pvalue_cwp_h2 = corr$p.value



############################################################
# CBP vs Healthy subgroup 
############################################################
if (!file.exists(file_cbp)){
  
  # All the healthy participants that were not part of the Tiemann study (Tiemann et. al, 2012) share similar characteristics in the recordings
  # original ID is "sub-FMpaXX"
  hc_nottiemann <- hc[!grepl("sub-FMhc[0-9]{2}",hc$original_id),]
  
  # Select a subset of healthy participants that match the characteristics of the CWP participants
  # i.e. that belonged to the same initial project and that don't differ in age or gender
  # We follow an iterative procedure that selects a random subset of healthy participants 
  # until we get evidence that there's no difference between groups in age and gender (BF<0.33)
  
  bf_age <- Inf
  bf_gender <- Inf
  selected_controls <- NULL
  counter <- 0
  while(bf_age >= 0.33 | bf_gender >= 0.33 & counter < 100){
    
    hc_cbp <- hc_nottiemann[sample(nrow(hc_nottiemann), nrow(cbp), replace = FALSE),]
    
    bf_results <- compute_bf(cbp, hc_cbp)
    
    bf_age <- bf_results$bf_age
    bf_gender <- bf_results$bf_gender
    counter <- counter +1
  }
  print("Exit Loop")
  # Save data in a csv
  data_cbp <- rbind(hc_cbp,cbp)
  
  # Regress out age from the aperiodic exponents in the subgroup of participants
  model <- lm(data_cbp$exp_PFC ~ data_cbp$age)
  res_apexp <- resid(model)
  res_cbp <-res_apexp[data_cbp$group == "pa"]
  res_hc_cbp <-res_apexp[data_cbp$group == "hc"]
  data_cbp$res_apexp_h1 = res_apexp
  
  # Regress out age from the aperiodic exponents in the cbp patients
  cbp <- data_cbp[data_cbp$group == "pa",]
  model <- lm(cbp$exp_PFC ~ cbp$age)
  res_apexp_cbp <-resid(model)
  # Regress out age from pain ratings in the cwp patients
  model <- lm(cbp$avg_pain_new ~ cbp$age)
  res_pain_cbp <-resid(model)
  data_cbp$res_apexp_h2[data_cbp$group == "pa"] = res_apexp_cbp
  data_cbp$res_pain_h2[data_cbp$group == "pa"] = res_pain_cbp
  
  write.table(data_cbp, file = "/rechenmagd3/Experiments/2023_1overf/results/statistics/e4_cbp.csv", sep=",", row.names = FALSE)
  
}else{
  data_cbp <- read.csv(file_cbp, header = TRUE, sep = ",")
  # subset_hc <- data_cbp[data_cbp$group == "hc",]
  # cbp <- data_cbp[data_cbp$group == "pa",]
}


# H1. Bayesian two sided independent samples t-test between groups on the residuals
bttest = ttestBF(res_hc_cbp, res_cbp, mu = 0, paired = FALSE, rscale = "medium", posterior = FALSE)
post = posterior(bttest,iterations = 1000)
bf_cbp_h1 = as.data.frame(bttest)$bf
postdelta_cbp_h1 = median(post[,"delta"])
# Frequentist t-test and effect size (Cohen's d)
fttest = t.test(res_hc_cbp, res_cbp, alternative = "two.sided", var.equal = FALSE)
pvalue_cbp_h1 = fttest$p.value
eff = cohen.d(res_hc_cbp,res_cbp)
cohens_d_cbp_h1 = eff$estimate

# H2. Bayesian correlation between age-corrected aperiodic exponents and age-corrected pain ratings
# in the CWP dataset
bcorrelation = correlationBF(y = res_apexp_cbp, x = res_pain_cbp, rscale = "ultrawide")
post = posterior(bcorrelation,iterations = 1000)
bf_cbp_h2 = as.data.frame(bcorrelation)$bf
postrho_cbp_h2 = median(post[,"rho"])
# Pearson's correlation
corr = cor.test(y = res_apexp_cbp, x = res_pain_cbp)
r_cbp_h2 = as.data.frame(corr$estimate)$cor
pvalue_cbp_h2 = corr$p.value





