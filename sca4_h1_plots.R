# Plot specification curve
#
# Cristina Gil, Elisabeth May, TUM, 16.02.2024

rm(list=ls()) 

# load libraries
library(specr)
library(tidyverse)
library(dplyr)
library(ggplot2)
library(cowplot)
library(viridis)

##################################

# load in data and select results for current electrode selection
filename = "C:\\Users\\Cristina\\1overf\\results\\sca\\specifications.txt"
specifications = read.csv(filename, header = TRUE, sep = ",")
filename = "C:\\Users\\Cristina\\1overf\\results\\sca\\randomizations_h1\\stats_orig.csv"
stats = read.csv(filename, header = TRUE, sep = ",")
results <- cbind(specifications,stats)

# load data from inference
filename <- "C:\\Users\\Cristina\\1overf\\results\\sca\\sca_h1_inference.csv"
inference <- read.csv(filename, header = TRUE, sep = ",")

########## prepare data ############################

# sort results based on effect size
results <- results %>%
  dplyr::arrange((.data$d))

# add specification number
results <- results %>%
  dplyr::mutate(specifications = 1:n())

# Data from the highlighted specification (main analysis)
highlightSpec = 1;
spec_highlighted = which(results$spec_id %in% c(1))
d_highlighted = results$d[spec_highlighted]
bf_highlighted = results$BF[spec_highlighted]

####################################################


########### plot curve   ###########################
p1 <-  results %>%
  ggplot(aes(x = specifications, y = d, ymin = d, ymax = d)) +
  geom_hline(yintercept= 0, linetype = "solid", color = "black", linewidth = 0.3) +
  geom_point(aes(color = log10(BF)), size = 0.75) +
  geom_point(aes(x = spec_highlighted, y = d_highlighted, color = log10(bf_highlighted)), size = 4, shape = 124) +
  scale_color_viridis(limits = c(-1, 1),option = "viridis", breaks = c(-1, log10(1/3), 0, log10(3),1)) +
  theme_minimal() +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_line(linewidth = 0.25)) + 
  labs(x = "", y = "Cohen's d") +
  ylim(-0.2, 0.2)  

#####################################################


########### plot choices ###########################
choices = c("epoch_length", "taper", "average_psd", "fooof_range", "fooof_knee")

value <- key <- NULL
p2 <- results %>%
  tidyr::gather(key, value, choices) %>% 
  dplyr::mutate(key = factor(.data$key, levels = choices)) %>% 
  ggplot(aes(x = .data$specifications, y = .data$value)) + 
  geom_point(aes(x = .data$specifications, y = .data$value, color = log10(.data$BF)), shape = 124, size = 3.35) + 
  scale_color_viridis(limits = c(-1, 1), option = "viridis") +
  theme_minimal() + 
  facet_grid(.data$key ~ 1, scales = "free_y", space = "free_y") + 
  theme(axis.line.y = element_line("black", size = 0.5), legend.position = "none", 
        panel.spacing = unit(0.75, "lines"), axis.text = element_text(colour = "black", size = 10), 
        panel.grid.major = element_line(size = 0.25),
        strip.text.x = element_blank(),
        strip.text.y = element_text(angle = 360, hjust = 0),
        panel.grid.minor = element_blank()) + 
  labs(x = "", y = "")

#####################################################

########### plot inference ##########################
p3 <- ggplot() +
  geom_line(aes(x = 1:48, y= inference$d_orig), color = "#440154") +
  geom_line(aes(x = 1:48, y= inference$d_median), color = "#808080") +
  geom_line(aes(x = 1:48, y= inference$lowPC), color = "#808080", linetype = "dashed") +
  geom_line(aes(x = 1:48, y= inference$highPC), color = "#808080", linetype = "dashed") +
  geom_hline(yintercept= 0, linetype = "solid", color = "black", size = 0.3) +
  theme_minimal() +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_line(linewidth = 0.25)) +
  labs(x = "specifications (nr), sorted by effect size", y = "Cohen's d") +
  ylim(-0.4, 0.4)

########### join plots and save as pdf #####################
pdf("C:\\Users\\Cristina\\1overf\\results\\sca\\sca_h1_plots.pdf")
cowplot::plot_grid(p1, p2, p3, align = "v",
                   axis = "rbl", rel_heights = c(2, 3.5, 1.5), ncol = 1)
dev.off()  # Close the PDF device
