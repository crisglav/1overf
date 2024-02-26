# Randomizations for the specification curve analysis of exploratory analysis 1.
#
# We randomize the age 500 times, so that age does not match the participants.
# For each randomization and specification, we perform a bayesian correlation 
# between age and aperiodic exponents.

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
age_original = participants$age

# Folder with the output data
scapath <- "/rechenmagd3/Experiments/2023_1overf/results/sca"
outpath <- file.path(scapath,"randomizations_e1")
if(!dir.exists(outpath)){
  dir.create(outpath)
}

# Load original specification curve 
mat_data <- readMat(file.path(scapath,"specs_ap_exp.mat"))

# Number of specifications, and participants
nSpec = dim(mat_data$exp)[1]
nSubj =  dim(mat_data$exp)[2]

# hardcoded number of randomizations 
nRand = 500

for (iRand in 0:nRand) {
  
  if (iRand == 0){
    # Original specification curve
    age = age_original
  }else{
    # Randomize age
    random_ix = sample(1:nSubj,nSubj)
    age <- age_original[random_ix]
  }
  
  # Initialize vectors to store results
  BF <- numeric(nSpec)
  PostRho <- numeric(nSpec)
  R <- numeric(nSpec)
  pvalue <- numeric(nSpec)
  
  for (iSpec in 1:nSpec) {
    # Aperiodic exponents for this spcification
    apexp = mat_data$exp[iSpec,]
    # Handle specifications that contain nan values for aperiodic exponents at certain participants
    participants_mask = !is.na(apexp)
    
    # Bayesian correlation between the ap. exp. and age
    bf = correlationBF(y = apexp[participants_mask], x = age[participants_mask], rscale = "ultrawide")
    post = posterior(bf,iterations = 1000)
    BF[iSpec] = as.data.frame(bf)$bf
    PostRho[iSpec] = median(post[,"rho"])
    
    # Pearson's correlation
    corr = cor.test(y = apexp[participants_mask], x = age[participants_mask])
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
