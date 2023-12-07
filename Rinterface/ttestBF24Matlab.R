# clear variables
#rm(list=ls())


# load libraries
#install.packages("BayesFactor")
library(BayesFactor)
# set working directory
setwd('/rechenmagd3/Experiments/2023_1overf/results/sca/Rinterface/')

partTab <- read.csv(file = 'tempTTest.csv', head=T)
datX = partTab$X[!is.na(partTab$X)]
datY = partTab$Y[!is.na(partTab$Y)]

TTestBF = ttestBF(datX, datY, mu = 0, paired = FALSE, rscale = "medium", posterior = FALSE)
TTestBFdf = as.data.frame(TTestBF)
BF = TTestBFdf$bf

TTestBFSummary = summary(ttestBF(datX, datY, mu = 0, paired = FALSE, rscale = "medium", posterior = TRUE, iteration = 1000))

PostDeltaXminY = TTestBFSummary$statistics[4,1]
resTab = data.frame(BF, PostDeltaXminY)

write.csv(resTab, 'tempTTestRes.csv')