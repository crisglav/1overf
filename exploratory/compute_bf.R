library(BayesFactor)

# load in data and select results for current electrode selection
filename = "C:\\Users\\Mitarbeiter\\1overf\\results\\statistics\\exp_PFC_real.csv"
results = read.csv(filename, header = TRUE, sep = ",")

res = results[,c("exp_PFC","group")]
boxplot(exp_PFC ~ group, data = res, main = "Aperiodic exponents")

# frequentist independent samples ttest
t.test(exp_PFC ~ group, data = res, var.eq=TRUE)

# bayesian independent samples ttest
bf = ttestBF(formula = exp_PFC ~ group, data = res)
bf