# Randomizations for the specification curve analysis of hypothesis 1.
#
# We randomize group labels 500 times. 
# For each randomization and specification, we perform a bayesian two-sided t-test 
# between groups, after correcting the aperiodic components for age.


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

# Load age data from the main analysis' results.
file <- "/rechenmagd3/Experiments/2023_1overf/results/statistics/exp_PFC_real.csv"
participants = read.csv(file, header = TRUE, sep = ",")
age = participants$age

# Folder with the output data
scapath <- "/rechenmagd3/Experiments/2023_1overf/results/sca"
outpath <- file.path(scapath,"randomizations_h1")
if(!dir.exists(outpath)){
  dir.create(outpath)
}

# Load aperiodic exponents for all the specifications
mat_data <- readMat(file.path(scapath,"specs_ap_exp.mat"))

# Patient mask
pa_mask_original = as.logical(mat_data$pa.mask)

# Number of specifications, and participants
nSpec = dim(mat_data$exp)[1]
nSubj =  dim(mat_data$exp)[2]
  
# hardcoded number of randomizations 
nRand = 500

for (iRand in 0:nRand) {
  
  if (iRand == 0){
    # Original specification curve
    pa_mask <- pa_mask_original
  }else{
    # Create a random mask for patients and healthy (patients are ones)
    pa_mask <- sample(rep(c(TRUE,FALSE), c(sum(pa_mask), sum(!pa_mask))))
  }
  
  # Initialize vectors to store results
  BF <- numeric(nSpec)
  PostDelta <- numeric(nSpec)
  d <- numeric(nSpec)
  pvalue <- numeric(nSpec)

  for (iSpec in 1:nSpec) {
    # Aperiodic exponents
    apexp = mat_data$exp[iSpec,]
    # Handle specifications that contain nan values for aperiodic exponents at certain participants
    participants_mask = !is.na(apexp)

    # --- OPTION 3. Regress out age from the apexp and then do a bayesian ttest on the residuals.
    model <- lm(apexp[participants_mask] ~ age[participants_mask])
    residuals_apexp <- rep(NA,nSubj)
    residuals_apexp[participants_mask] <- resid(model)

    res_pa = residuals_apexp[pa_mask & participants_mask]
    res_hc  = residuals_apexp[!pa_mask & participants_mask]

    # Bayesian t-test between ap. exp. residuals of patients and healthy participans
    bf = ttestBF(res_hc, res_pa, mu = 0, paired = FALSE, rscale = "medium", posterior = FALSE)
    post = posterior(bf,iterations = 1000)
    BF[iSpec] = as.data.frame(bf)$bf
    PostDelta[iSpec] = median(post[,"delta"])

    # Frequentist t-test
    ttest = t.test(res_hc,res_pa, alternative = "two.sided", var.equal = FALSE)
    pvalue[iSpec] = ttest$p.value

    # Effect size
    eff = cohen.d(res_hc,res_pa)
    d[iSpec] = eff$estimate
    
    # # --- OPTION 2. Linear regression. 
    # # Problem: which BF for the model I use. The ones that are output are relative to each other. How do I sample from the group model to get a posterior delta>?
    # # ANCOVA, operationalized as a linear regression with group as a dummy variable (pa_mask), and age as a regressor
    # df <- data.frame(apexp, pa_mask, age)
    # bf <- regressionBF(apexp ~ pa_mask + age, data <- df)
    
    # # --- OPTION 1. Bayesian ttest, I don't take into account age
    # pa_exp <- apexp[pa_mask & participants_mask]
    # hc_exp <- apexp[!pa_mask & participants_mask]
    # 
    # # Bayesian t-test
    # bf = ttestBF(hc_exp, pa_exp, mu = 0, paired = FALSE, rscale = "medium", posterior = FALSE)
    # post = posterior(bf, iterations = 1000)
    # BF[iSpec] = as.data.frame(bf)$bf
    # PostDelta[iSpec] = median(post[,"delta"])
    # 
    # # Frequentist t-test
    # ttest = t.test(hc_exp,pa_exp, alternative = "two.sided", var.equal = FALSE)
    # pvalue[iSpec] = ttest$p.value
    # 
    # # Effect size
    # eff = cohen.d(hc_exp,pa_exp)
    # d[iSpec] = eff$estimate
  }
  
  resTab = data.frame(BF,PostDelta,d,pvalue)
  
  # Export dataframe as csv
  if (iRand == 0){
    csvname <- "stats_orig.csv"
  }else{
    csvname <- sprintf("stats_rand%0.3d.csv",iRand)
  }
  csvpath <- file.path(outpath,csvname)
  write.csv(resTab, file = csvpath, row.names = FALSE)
}
