# Relationship betwen age and aperiodic exponents
#
# Cristina Gil, TUM, 26.02.2024

rm(list=ls())
library(BayesFactor)

# Load prepared data from csv
datapath <- "/rechenmagd3/Experiments/2023_1overf/results/statistics/"
# datapath <- "C:\\Users\\Cristina\\1overf\\results\\statistics"
file <- paste("exp_PFC_real", ".csv", sep = "")
participants = read.csv(file.path(datapath,file), header = TRUE, sep = ",")

#################################################
# Exploratory 1: Are aperiodic exponents in the mPFC correlated with age?
# Bayesian correlation between aperiodic exponents and age
bcorrelation = correlationBF(y = participants$exp_PFC, x = participants$age, rscale = "ultrawide")
post = posterior(bcorrelation,iterations = 1000)
bf = as.data.frame(bcorrelation)$bf
postrho = median(post[,"rho"])

# Pearson's correlation
corr = cor.test(y = participants$exp_PFC, x = participants$age)
R = as.data.frame(corr$estimate)$cor
