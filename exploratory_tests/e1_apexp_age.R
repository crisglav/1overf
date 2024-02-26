# Relationship betwen age and aperiodic exponents
#
# Cristina Gil, TUM, 26.02.2024

rm(list=ls())
library(BayesFactor)
library(ggplot2)
library(svglite)

# Load prepared data from csv
# datapath <- "/rechenmagd3/Experiments/2023_1overf/results/statistics/"
datapath <- "C:\\Users\\Cristina\\1overf\\results\\statistics"
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

# Plot
p1 <-
  ggplot(participants, aes(x = age, y = exp_PFC)) + 
  geom_point(color = "#A7C133") + 
  geom_smooth(method = "lm", se = TRUE, color = "grey", alpha = 0.2) +
  theme_classic() +
  labs(title = "Relationship between age and aperiodic, all participants",
       x = "Age",
       y = "Aperiodic exponent")

# Save the plot as an SVG file
fig_path <- "C:\\Users\\Cristina\\1overf\\results\\figures\\e1_apexp_age.svg"
ggsave(fig_path, p1, device = "svg", width = 6, height = 4)