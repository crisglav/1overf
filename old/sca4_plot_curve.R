# Plot specification curve
#
# Cristina Gil, Elisabeth May, TUM, 17.01.2024

rm(list=ls()) 
setwd("C:/Users/Mitarbeiter/Downloads")

# load libraries
library(specr)
library(tidyverse)
library(dplyr)
library(ggplot2)
library(cowplot)

##################################

# load in data and select results for current electrode selection
filename = "C:\\Users\\Mitarbeiter\\1overf\\results\\sca\\specifications.txt"
specifications = read.csv(filename, header = TRUE, sep = ",")
filename = "C:\\Users\\Mitarbeiter\\1overf\\results\\sca\\Rinterface\\stats_orig.csv"
stats = read.csv(filename, header = TRUE, sep = ",")
results <- cbind(specifications,stats)

########## prepare data ############################

# sort results based on effect size, add specification number and define colors according to BF10
results <- results %>% 
  dplyr::arrange((.data$d))  
results <- results %>% 
  dplyr::mutate(specifications = 1:n(), 
                color = case_when(BF > 3 ~ "#377eb8", 
                                  BF < .3333 ~ "#e41a1c",
                                  TRUE ~ "darkgrey"))

# select data of  highlighted specifications for highlighting in plot
highlightSpec = 1;
indexhighlightedSpecs = which(results$spec_id %in% c(1))
d = results$d[indexhighlightedSpecs]
# confLow = results$ci_inf[indexhighlightedSpecs] # fake confidence intervals
# confHigh = results$ci_sup[indexhighlightedSpecs] # fake confidence intervals
highlightedData <- data.frame(d, indexhighlightedSpecs)

####################################################


########### plot curve   ###########################
p1 <-  results %>%
  ggplot(aes(x = .data$specifications, y = .data$d, 
             ymin = .data$d, ymax = .data$d, color = .data$color)) +
  geom_hline(yintercept= 0, linetype = "solid", color = "black", size = 0.3) +
  geom_point(aes(color = .data$color), size = 0.4) +
  theme_minimal() + 
  scale_color_identity() + 
  theme(strip.text = element_blank(), axis.line = element_line("black", size = 0.5),  # use linewidth instead of size for higher ggplot2 versions
        legend.position = "none", 
        panel.spacing = unit(0.75, "lines"), 
        axis.text = element_text(colour = "black")) + 
  #geom_pointrange(alpha = 0.5, size = 0.1, fatten = 0.1) + 
  # mark highlighted data
  geom_point(data = highlightedData, 
             aes(x = .data$indexhighlightedSpecs, y = .data$d),
             color = "black", size = 1.5, shape = 16) +
  labs(x = " ", y = "Cohen's d")+
  ylim(-0.3, 0.3)  

#####################################################


########### plot choices ###########################
choices = c("epoch_length", "taper", "average_psd", "fooof_range", "fooof_knee")

aux <- results %>%
  tidyr::gather(key, value, choices) %>% 
  dplyr::mutate(key = factor(.data$key, levels = choices))

value <- key <- NULL
p2 <- results %>%
  tidyr::gather(key, value, choices) %>% 
  dplyr::mutate(key = factor(.data$key, levels = choices)) %>% 
  ggplot(aes(x = .data$specifications, y = .data$value, color = .data$color)) + 
  geom_point(aes(x = .data$specifications, y = .data$value), shape = 124, size = 3.35) + 
  scale_color_identity() + 
  theme_minimal() + 
  facet_grid(.data$key ~ 1, scales = "free_y", space = "free_y") + 
  theme(axis.line = element_line("black", size = 0.5), legend.position = "none", 
        panel.spacing = unit(0.75, "lines"), axis.text = element_text(colour = "black", size = 10), 
        strip.text.x = element_blank(),
        strip.text.y = element_text(angle = 360, hjust = 0)) + 
  labs(x = "", y = "")

#####################################################


########### join plots and save as pdf #####################
pdf("C:\\Users\\Mitarbeiter\\1overf\\results\\sca\\results_1overf_specr_curve.pdf")
cowplot::plot_grid(p1, p2, align = "v",
                   axis = "rbl", rel_heights = c(2, 3.5, 1.5), ncol = 1)
dev.off()  # Close the PDF device
