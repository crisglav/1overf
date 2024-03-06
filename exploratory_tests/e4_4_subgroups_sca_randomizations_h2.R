# Subgroups test
# Randomizations for the specification curve analysis of hypothesis 2 (correlation
# between pain intensity and aperiodic exponents in patients).
#
# We randomize the order of pain ratings 500 times, so that pain ratings do not match the participants. 
# For each randomization and specification, we perform a bayesian correlation 
# between age-corrected aperiodic exponents and age-corrected randomized pain ratings.
# Aperiodic exponents and pain ratings are corrected for age by regressing out age from them.


# clear variables, set seed
rm(list=ls())
set.seed(123)

# load libraries
# install.packages("BayesFactor")
# install.packages("R.matlab")
# install.packages("effsize")
library(BayesFactor)
library(R.matlab)
library(effsize)

# Load aperiodic exponents for all the specifications
scapath <- "/rechenmagd3/Experiments/2023_1overf/results/sca/"
mat_data <- readMat(file.path(scapath,"specs_ap_exp.mat"))

# Create output folder
outpath <- file.path(scapath,"randomizations_e4_h2")
if(!dir.exists(outpath)){
  dir.create(outpath)
}
# Load participants file with all the participants. We need this file to create a mask, so
# that we can extract from the aperiodic exponents for all the specifications the subgroups
# of participants
file_participants <- "/rechenmagd3/Experiments/2023_1overf/results/statistics/exp_PFC_real.csv"
participants = read.csv(file_participants, header = TRUE, sep = ",")

# Load data from the subgroups of patients and healthy participants.
file_cwp <- "/rechenmagd3/Experiments/2023_1overf/results/statistics/e4_cwp.csv"
data_cwp <- read.csv(file_cwp, header = TRUE, sep = ",")
file_cbp <- "/rechenmagd3/Experiments/2023_1overf/results/statistics/e4_cbp.csv"
data_cbp <- read.csv(file_cbp, header = TRUE, sep = ",")

# Mask of all the participants of the subgroup (cwp and matched healthy to the cwp group)
mask_cwp <- participants$participant_id %in% data_cwp$participant_id & participants$group == "pa"
mask_cbp <- participants$participant_id %in% data_cbp$participant_id & participants$group == "pa"

# Extract aperiodic exponents, pain ratings and age for the subgroup
apexp_cwp <- mat_data$exp[,mask_cwp]
age_cwp <- participants$age[mask_cwp]
pain_cwp <- data_cwp$avg_pain_new[data_cwp$group == "pa"]

apexp_cbp <- mat_data$exp[,mask_cbp]
age_cbp <- participants$age[mask_cbp]
pain_cbp <- data_cbp$avg_pain_new[data_cbp$group == "pa"]

# Regress out age from pain ratings
model <- lm(pain_cwp ~ age_cwp)
pain_res_cwp_original <- resid(model)
model <- lm(pain_cbp ~ age_cbp)
pain_res_cbp_original <- resid(model)

# Number of specifications, participants of the subgroup and randomizations
nSpec = dim(apexp_cwp)[1]
nSubj_cwp =  dim(apexp_cwp)[2]
nSubj_cbp =  dim(apexp_cbp)[2]
nRand = 500

for (iRand in 0:nRand) {
  # Initialize dataframe to store results
  results = matrix(nrow = nSpec, ncol = 8)
  colnames(results) <- c("BF_cwp","PostR_cwp","pvalue_cwp","R_cwp","BF_cbp","PostR_cbp","pvalue_cbp","R_cbp")
  
  if (iRand == 0){
    # Original mask
    pain_res_cwp <- pain_res_cwp_original
    pain_res_cbp <- pain_res_cbp_original
  }else{
    # Randomize pain ratings so that they are no longer related to the original participant
    random_ix = sample(1:nSubj_cwp,nSubj_cwp)
    pain_res_cwp <- pain_res_cwp_original[random_ix]
    random_ix = sample(1:nSubj_cbp,nSubj_cbp)
    pain_res_cbp <- pain_res_cbp_original[random_ix]
  }
  
  for (iSpec in 1:nSpec) {
    
    ########## CWP ###########
    # Aperiodic exponents of patients
    apexp = apexp_cwp[iSpec,]
    # Handle specifications that contain nan values for aperiodic exponents at certain participants
    participants_mask = !is.na(apexp)
    
    # Regress out age from apexp
    model <- lm(apexp[participants_mask] ~ age_cwp[participants_mask])
    residuals_apexp <- resid(model)
    
    # Bayesian correlation between the ap. exp. residuals and pain ratings residuals
    bf = correlationBF(y = residuals_apexp, x = pain_res_cwp[participants_mask], rscale = "ultrawide")
    post = posterior(bf,iterations = 1000)
    results[iSpec,"BF_cwp"] = as.data.frame(bf)$bf
    results[iSpec,"PostR_cwp"] = median(post[,"rho"])
    
    # Pearson's correlation
    corr = cor.test(y = residuals_apexp, x = pain_res_cwp[participants_mask])
    results[iSpec,"R_cwp"] = as.data.frame(corr$estimate)$cor
    results[iSpec,"pvalue_cwp"] = corr$p.value
    
    ########## CBP ###########
    # Aperiodic exponents of patients
    apexp = apexp_cbp[iSpec,]
    # Handle specifications that contain nan values for aperiodic exponents at certain participants
    participants_mask = !is.na(apexp)
    
    # Regress out age from apexp
    model <- lm(apexp[participants_mask] ~ age_cbp[participants_mask])
    residuals_apexp <- resid(model)
    
    # Bayesian correlation between the ap. exp. residuals and pain ratings residuals
    bf = correlationBF(y = residuals_apexp, x = pain_res_cbp[participants_mask], rscale = "ultrawide")
    post = posterior(bf,iterations = 1000)
    results[iSpec,"BF_cbp"] = as.data.frame(bf)$bf
    results[iSpec,"PostR_cbp"] = median(post[,"rho"])
    
    # Pearson's correlation
    corr = cor.test(y = residuals_apexp, x = pain_res_cbp[participants_mask])
    results[iSpec,"R_cbp"] = as.data.frame(corr$estimate)$cor
    results[iSpec,"pvalue_cbp"] = corr$p.value
  }
  
  resTab = data.frame(results)
  
  # Export dataframe as csv
  if (iRand == 0){
    csvname <- "stats_orig.csv"
  }else{
    csvname <- sprintf("stats_rand%0.3d.csv",iRand)
  }
  csvpath <- file.path(outpath,csvname)
  write.csv(resTab, file = csvpath, row.names = FALSE)
}
