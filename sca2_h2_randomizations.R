# Randomizations for the specification curve analysis of hypothesis 2.
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

# Load age and pain ratings. File already in the stats folder because it was adjusted to fill missing values of avg_pain rating with current pain ratings
file <- "/rechenmagd3/Experiments/2023_1overf/results/statistics/patients_PFC_real.csv"
patients = read.csv(file, header = TRUE, sep = ",")
pain_mask = !is.na(patients$avg_pain)
pain_ratings_original = patients$avg_pain[pain_mask]
age = patients$age[pain_mask]

# Regress out age from pain ratings
model_pain <- lm(pain_ratings_original ~ age)
residuals_pain_original <- resid(model_pain)

# Folder with the output data
scapath <- "/rechenmagd3/Experiments/2023_1overf/results/sca"
outpath <- file.path(scapath,"randomizations_h2")
if(!dir.exists(outpath)){
  dir.create(outpath)
}

# Load original specification curve 
mat_data <- readMat(file.path(scapath,"specs_ap_exp.mat"))

# Patient mask
pa_mask = as.logical(mat_data$pa.mask)

# Extract the aperiodic exponents of patients that have pain ratings
apexp_pa = mat_data$exp[,pa_mask]
apexp_pa = apexp_pa[,pain_mask]

# Number of specifications, and participants
nSpec = dim(apexp_pa)[1]
nSubj =  dim(apexp_pa)[2]

# hardcoded number of randomizations 
nRand = 500

for (iRand in 0:nRand) {
  
  if (iRand == 0){
    # Original specification curve
    residuals_pain = residuals_pain_original
  }else{
    # Randomize pain ratings so that they are no longer related to the original participant
    random_ix = sample(1:nSubj,nSubj)
    residuals_pain <- residuals_pain_original[random_ix]
  }
  
  # Initialize vectors to store results
  BF <- numeric(nSpec)
  PostRho <- numeric(nSpec)
  R <- numeric(nSpec)
  pvalue <- numeric(nSpec)
  
  for (iSpec in 1:nSpec) {
    # Aperiodic exponents of patients
    apexp = apexp_pa[iSpec,]
    # Handle specifications that contain nan values for aperiodic exponents at certain participants
    participants_mask = !is.na(apexp)

    # Regress out age from apexp
    model_apexp <- lm(apexp[participants_mask] ~ age[participants_mask])
    residuals_apexp <- resid(model_apexp)
    
    # Bayesian correlation between the ap. exp. residuals and pain ratings residuals
    bf = correlationBF(y = residuals_apexp, x = residuals_pain[participants_mask], rscale = "ultrawide")
    post = posterior(bf,iterations = 1000)
    BF[iSpec] = as.data.frame(bf)$bf
    PostRho[iSpec] = median(post[,"rho"])
    
    # Pearson's correlation
    corr = cor.test(y = residuals_apexp, x = residuals_pain[participants_mask])
    R[iSpec] = as.data.frame(corr$estimate)$cor
    pvalue[iSpec] = corr$p.value

  }
  
  resTab = data.frame(BF,PostRho,R,pvalue)
  
  # Export dataframe as csv
  if (iRand == 0){
    csvname <- "stats_orig.csv"
  }else{
    csvname <- sprintf("stats_rand%0.3d.csv",iRand)
  }
  csvpath <- file.path(outpath,csvname)
  write.csv(resTab, file = csvpath, row.names = FALSE)
}
