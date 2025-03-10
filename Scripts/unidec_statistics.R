
library(ggplot2)
library(dplyr)
library(tidyverse)
library(ggpubr)


data <- read.csv("C:/Users/ives435/OneDrive - PNNL/Desktop/GPCR paper/unidec_stats.csv")

palette <- c("darkgreen", "darkred", "darkblue")
boxsize <- 1
psize <- 8

#unmod 
unmod <- data %>%
   filter(Proteoform == 1)

my_comparisons <- list(
   c("1", "2"), c("2", "3"),
   c("1", "3")
)

unm <- ggboxplot(size = boxsize,
   unmod, x = "Peak", y = "Height",
   color = "Peak", palette = palette,
   add = "dotplot"
)+
   stat_compare_means(comparisons = my_comparisons, method = "t.test", label = "p.signif", size=psize)+
   labs(x="Migration time region", y = "Relative Intensity (UniDec)", title="Unmodified", color = "Region", fill="Region")+
   theme_classic(base_size = 24)+
   ylim(50,120)

#monophos
monophos <- data %>%
   filter(Proteoform == 2)

my_comparisons <- list(
   c("1", "2"), c("2", "3"),
   c("1", "3")
)

mp <- ggboxplot(size = boxsize,
   monophos, x = "Peak", y = "Height",
   color = "Peak", palette = palette,
   add = "dotplot"
)+
   stat_compare_means(comparisons = my_comparisons, method = "t.test", label = "p.signif", size=psize)+
   labs(x="Migration time region", y = "Relative Intensity (UniDec)", title="Monophosphorylated", color = "Region", fill="Region")+
   theme_classic(base_size = 24)+
   ylim(95,105)


#diphos
diphos <- data %>%
   filter(Proteoform == 3)

my_comparisons <- list(
   c("1", "2"), c("2", "3"),
   c("1", "3")
)

dp <- ggboxplot(size = boxsize,
   diphos, x = "Peak", y = "Height",
   color = "Peak", palette = palette,
   add = "dotplot"
)+
   stat_compare_means(comparisons = my_comparisons, method = "t.test", label = "p.signif", size=psize)+
   labs(x="Migration time region", y = "Relative Intensity (UniDec)", title="Diphosphorylated", color = "Region", fill="Region")+
   theme_classic(base_size = 24)+
   ylim(25,55)
dp


#diphos
acyl <- data %>%
   filter(Proteoform == 4)

my_comparisons <- list(
   c("1", "2"), c("2", "3"),
   c("1", "3")
)

ac <- ggboxplot(size = boxsize,
   acyl, x = "Peak", y = "Height",
   color = "Peak", palette = palette,
   add = "dotplot"
)+
   stat_compare_means(comparisons = my_comparisons, method = "t.test", label = "p.signif", size=psize)+
   labs(x="Migration time region", y = "Relative Intensity (UniDec)", title="Palmitoyl", color = "Region", fill="Region")+
   theme_classic(base_size = 24)+
   ylim(70,95)
ac

library(patchwork)

combined_plot <- (unm + mp) / (dp + ac) +
   plot_annotation(tag_levels = 'A')
combined_plot

ggsave(combined_plot, filename = "C:/Users/ives435/OneDrive - PNNL/Desktop/GPCR paper/unidecstats.png",
       scale = 2, width = 7, height = 7, units = "in", dpi = 600)


