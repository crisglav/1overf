# Subgroups test
# Randomizations for the specification curve analysis of hypothesis 1 (differences between 
# matched healthy and subgroup of patients).
#
# We randomize group labels 500 times. 
# For each randomization and specification, we perform a bayesian two-sided t-test 
# between groups, after correcting the aperiodic exponents for age.


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
outpath <- file.path(scapath,"randomizations_e4_h1")
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
mask_cwp_all <- participants$participant_id %in% data_cwp$participant_id
cwp_all <- participants[mask_cwp_all,] 
mask_cwp <- cwp_all$group == "pa"
mask_cwp_hc <- cwp_all$group == "hc"

mask_cbp_all <- participants$participant_id %in% data_cbp$participant_id
cbp_all <- participants[mask_cbp_all,] 
mask_cbp <- cbp_all$group == "pa"
mask_cbp_hc <- cbp_all$group == "hc"

# Extract aperiodic exponents and age for the subgroup
apexp_cwp_all <- mat_data$exp[,mask_cwp_all]
age_cwp_all <- participants$age[mask_cwp_all]

apexp_cbp_all <- mat_data$exp[,mask_cbp_all]
age_cbp_all <- participants$age[mask_cbp_all]

# Number of specifications, participants of the subgroup and randomizations
nSpec = dim(apexp_cwp_all)[1]
nSubj_cwp =  dim(apexp_cwp_all)[2]
nSubj_cbp =  dim(apexp_cbp_all)[2]
nRand = 500

for (iRand in 0:nRand) {
  # Initialize dataframe to store results
  results = matrix(nrow = nSpec, ncol = 8)
  colnames(results) <- c("BF_cwp","PostDelta_cwp","pvalue_cwp","d_cwp","BF_cbp","PostDelta_cbp","pvalue_cbp","d_cbp")
  
  if (iRand == 0){
    # Original mask
    pa_mask_cwp <- mask_cwp
    pa_mask_cbp <- mask_cbp
  }else{
    # Create a random mask for patients and healthy (patients are ones)
    pa_mask_cwp <- sample(rep(c(TRUE,FALSE), c(sum(mask_cwp), sum(mask_cwp_hc))))
    pa_mask_cbp <- sample(rep(c(TRUE,FALSE), c(sum(mask_cwp), sum(mask_cwp_hc))))
  }
  
  for (iSpec in 1:nSpec) {
    
    ########## CWP ###########
    # Aperiodic exponents
    apexp = apexp_cwp_all[iSpec,]
    # Handle specifications that contain nan values for aperiodic exponents at certain participants
    participants_mask = !is.na(apexp)
    
    # Regress out age from the apexp and then do a bayesian ttest on the residuals.
    model <- lm(apexp[participants_mask] ~ age_cwp_all[participants_mask])
    residuals_apexp <- rep(NA,nSubj_cwp)
    residuals_apexp[participants_mask] <- resid(model)
    
    res_pa = residuals_apexp[pa_mask_cwp & participants_mask]
    res_hc  = residuals_apexp[!pa_mask_cwp & participants_mask]
    
    # Bayesian t-test between ap. exp. residuals of patients and healthy participans
    bf = ttestBF(res_hc, res_pa, mu = 0, paired = FALSE, rscale = "medium", posterior = FALSE)
    post = posterior(bf,iterations = 1000)
    results[iSpec,"BF_cwp"] = as.data.frame(bf)$bf
    results[iSpec,"PostDelta_cwp"] = median(post[,"delta"])
    
    # Frequentist t-test
    ttest = t.test(res_hc,res_pa, alternative = "two.sided", var.equal = FALSE)
    results[iSpec,"pvalue_cwp"] = ttest$p.value
    
    # Effect size
    eff = cohen.d(res_hc,res_pa)
    results[iSpec,"d_cwp"] = eff$estimate
    
    ########## CBP ###########
    # Aperiodic exponents
    apexp = apexp_cbp_all[iSpec,]
    # Handle specifications that contain nan values for aperiodic exponents at certain participants
    participants_mask = !is.na(apexp)
    
    # Regress out age from the apexp and then do a bayesian ttest on the residuals.
    model <- lm(apexp[participants_mask] ~ age_cbp_all[participants_mask])
    residuals_apexp <- rep(NA,nSubj_cbp)
    residuals_apexp[participants_mask] <- resid(model)
    
    res_pa = residuals_apexp[pa_mask_cbp & participants_mask]
    res_hc  = residuals_apexp[!pa_mask_cbp & participants_mask]
    
    # Bayesian t-test between ap. exp. residuals of patients and healthy participans
    bf = ttestBF(res_hc, res_pa, mu = 0, paired = FALSE, rscale = "medium", posterior = FALSE)
    post = posterior(bf,iterations = 1000)
    results[iSpec,"BF_cbp"] = as.data.frame(bf)$bf
    results[iSpec,"PostDelta_cbp"] = median(post[,"delta"])
    
    # Frequentist t-test
    ttest = t.test(res_hc,res_pa, alternative = "two.sided", var.equal = FALSE)
    results[iSpec,"pvalue_cbp"] = ttest$p.value
    
    # Effect size
    eff = cohen.d(res_hc,res_pa)
    results[iSpec,"d_cbp"] = eff$estimate
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
