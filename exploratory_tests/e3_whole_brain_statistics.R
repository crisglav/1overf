# Statistics of the main research questions

rm(list=ls())
library(BayesFactor)
library(R.matlab)


# Load prepared data from matlab with aperiodic exponents and offsets and age
datapath <- "/rechenmagd3/Experiments/2023_1overf/results/features/fooof_matlab/whole_brain/e3_whole_brain.mat"
mat_data <- readMat(datapath)

# For each ROI regress out age from the aperiodic exponent and perform a independent sambples bayesian t-test
# between healthy and patients
nRoi = dim(mat_data$apexp)[2]
bf_apexp <- numeric(nRoi)
for (iRoi in 1:nRoi) {
  apexp = mat_data$apexp[,iRoi]
  
  # Regress out age from aperiodic exponents
  model_apexp <- lm(apexp ~ mat_data$age)
  res_apexp <- resid(model_apexp)
  
  # Separate residuals based on group
  res_pa = res_apexp[as.logical(mat_data$pa.mask)]
  res_hc = res_apexp[as.logical(mat_data$hc.mask)]
  
  # Bayesian two sided independent samples t-test between groups on the residuals
  bttest = ttestBF(res_hc, res_pa, mu = 0, paired = FALSE, rscale = "medium", posterior = FALSE)
  bf_apexp[iRoi] = as.data.frame(bttest)$bf
  
  # apexp_pa = apexp[mat_data$pa.mask]
  # apexp_hc = apexp[mat_data$hc.mask]
  
}

#################################################
# Hypothesis 1: Do aperiodic exponents in the mPFC differ between patients and healthy participants?

# Regress out age from aperiodic exponents
model <- lm(participants$exp_PFC ~ participants$age)
res_apexp <- resid(model)

# Separate residuals based on group
res_pa = res_apexp[as.logical(participants$group_binary)]
res_hc = res_apexp[!as.logical(participants$group_binary)]

# Bayesian two sided independent samples t-test between groups on the residuals
bttest = ttestBF(res_hc, res_pa, mu = 0, paired = FALSE, rscale = "medium", posterior = FALSE)
post = posterior(bttest,iterations = 1000)
bf_h1 = as.data.frame(bttest)$bf
postdelta_h1 = median(post[,"delta"])

# Frequentist t-test and effect size (Cohen's d)
fttest = t.test(res_hc,res_pa, alternative = "two.sided", var.equal = FALSE)
pvalue_h1 = fttest$p.value
eff = cohen.d(res_hc,res_pa)
cohens_d_h1 = eff$estimate

#################################################
# Hypothesis 2: Do aperiodic exponents in the mPFC correlate with pain ratings in patients?
# Discard the data from the patient that does not have a pain rating
pain_mask = !is.na(patients$avg_pain)
avg_pain = patients$avg_pain[pain_mask]
age_pa = patients$age[pain_mask]
apexp_pa = patients$exp_PFC[pain_mask]

# Regress out age from pain ratings and aperiodic exponents
model_pain <- lm(avg_pain ~ age_pa)
res_pain <- resid(model_pain)
model_apexp <- lm(apexp_pa ~ age_pa)
res_apexp_pa <- resid(model_apexp)

# Bayesian correlation between the residuals
bcorrelation = correlationBF(y = res_apexp_pa, x = res_pain, rscale = "ultrawide")
post = posterior(bcorrelation,iterations = 1000)
bf_h2 = as.data.frame(bcorrelation)$bf
postrho_h2 = median(post[,"rho"])

# Pearson's correlation
corr = cor.test(y = res_apexp_pa, x = res_pain)
r_h2 = as.data.frame(corr$estimate)$cor
pvalue_h2 = corr$p.value