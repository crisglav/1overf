# clear variables
rm(list=ls())

# load libraries
# install.packages("BayesFactor")
# install.packages("R.matlab")
# install.packages("effsize")
library(BayesFactor)
library(R.matlab)
library(effsize)

filepath <- "/rechenmagd3/Experiments/2023_1overf/results/sca"

# hardcoded number of randomizations and specifications
nRand = 500
nSpec = 48

for (iRand in 1:nRand) {
  
  if (iRand == 0){
    filename <- "specs_ap_exp.mat"
  }else{
    filename <- sprintf("specs_ap_exp_rand%0.3d.mat",iRand)
  }
  
    
  file <- file.path(filepath,filename)
  mat_data <- readMat(file)
  
  # Initialize vectors to store results
  BF <- numeric(nSpec)
  PostDelta <- numeric(nSpec)
  pvalue <- numeric(nSpec)
  tvalue <- numeric(nSpec)
  d <- numeric(nSpec)
  
  for (iSpec in 1:nSpec) {
    hc <- mat_data$exp[iSpec, as.logical(mat_data$hc.mask)]
    pa <- mat_data$exp[iSpec, as.logical(mat_data$pa.mask)]
    
    # remove nans
    hc <- hc[is.finite(hc)]
    pa <- pa[is.finite(pa)]
    
    # Linear regression, group is a dummy variable, age is a covariate
    
    # Bayesian t-test
    TTestBF = ttestBF(hc, pa, mu = 0, paired = FALSE, rscale = "medium", posterior = FALSE)
    TTestBFSummary = summary(ttestBF(hc, pa, mu = 0, paired = FALSE, rscale = "medium", posterior = TRUE, iteration = 1000))
    
    BF[iSpec] = as.data.frame(TTestBF)$bf
    PostDelta[iSpec] = TTestBFSummary$statistics[4,1]
    
    # Frequentist t-test
    ttest = t.test(hc,pa, alternative = "two.sided", var.equal = FALSE)
    pvalue[iSpec] = ttest$p.value
    tvalue[iSpec] = ttest$statistic
    
    # Effect size
    eff = cohen.d(hc,pa)
    d[iSpec] = eff$estimate
  }
  
  resTab = data.frame(BF, PostDelta,pvalue,tvalue,d)
  
  # Export dataframe as csv
  if (iRand == 0){
    csvname <- "stats_orig.csv"
  }else{
    csvname <- sprintf("stats_rand%0.3d.csv",iRand)
  }
  csvpath <- file.path("/rechenmagd3/Experiments/2023_1overf/results/sca/Rinterface",csvname)
  write.csv(resTab, file = csvpath, row.names = FALSE)
}
