# Statistics of the main research questions

rm(list=ls())
library(BayesFactor)
library(effsize)

# Analysis type (real or blinded)
analysis <- "real"

# Load prepared data from csv
datapath <- "/rechenmagd3/Experiments/2023_1overf/results/statistics/"
file_h1 <- paste("exp_PFC_", analysis, ".csv", sep = "")
file_h2 <- paste("patients_PFC_", analysis, ".csv", sep = "")
participants = read.csv(file.path(datapath,file_h1), header = TRUE, sep = ",")
patients = read.csv(file.path(datapath,file_h2), header = TRUE, sep = ",")

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