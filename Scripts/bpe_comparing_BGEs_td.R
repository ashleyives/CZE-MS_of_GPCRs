library(ggplot2)
library(zoo)
library(dplyr)
library(rawrr)

# File paths
file_paths <- c(
   "C:/Users/ives435/OneDrive - PNNL/Desktop/GPCR paper/massive/20220729_ani0287_cap52_newb2ar_pngasef30min37C_01_SID100.raw",
   "C:/Users/ives435/OneDrive - PNNL/Desktop/GPCR paper/massive/to_add/20220216_ani0287_b2ar_3hrPNGaseFRT_noZEBA_30kV_100secInj_MS1_04_3pAAnosulfolane.raw",
   "C:/Users/ives435/OneDrive - PNNL/Desktop/GPCR paper/massive/20220729_ani0287_cap52_newb2ar_pngasef30min37C_01_SID50_CID20_CS27_1850mz.raw"
)

# Initialize an empty data frame to store results
combined_data <- data.frame()

# Process each file
for (file in file_paths) {
   # Read chromatogram data for each file
   C <- rawrr::readChromatogram(file, type = "bpc")
   
   # Extract times and intensities
   times <- as.numeric(C$times)
   intensities <- as.numeric(C$intensities)
   
   # Create data frame and add file column
   data <- data.frame(times = times, intensities = intensities, file = basename(file))
   
   # Append to combined data
   combined_data <- rbind(combined_data, data)
}

# Check the levels of your 'file' factor
levels(combined_data$file)

# Create a named vector for your labels
facet_labels <- c(`20220216_ani0287_b2ar_3hrPNGaseFRT_noZEBA_30kV_100secInj_MS1_04_3pAAnosulfolane.raw` = "A - 3% Acetic Acid, SID = 100", 
                  `20220729_ani0287_cap52_newb2ar_pngasef30min37C_01_SID100.raw` = "B - 3% Acetic Acid+10%sulfolane, SID = 100", 
                  `20220729_ani0287_cap52_newb2ar_pngasef30min37C_01_SID50_CID20_CS27_1850mz.raw` = "C - 3% Acetic Acid+10%sulfolane, SID = 50") 

# Create a named vector for legend labels with empty strings to hide them
legend_labels <- c("groupA" = "")

# Plot
plot <- ggplot(combined_data, aes(x = times, y = intensities)) +
   geom_line(size = 1) +
   labs(y = "Rel. Abundance", x = "Migration time (min)") + # Change legend title here
   theme_classic(base_size = 24) +
   xlim(10, 35) +
   facet_wrap(~file, ncol = 1, labeller = as_labeller(facet_labels), scales="free") +
   theme(
      text = element_text(family = "Arial"), # Set font to Arial
      strip.text = element_text(face = "bold", family = "Arial", size = rel(1.2)), # Make facet label text bold and bigger
      strip.text.y.left = element_text(angle = 0, hjust = 0, face = "bold", family = "Arial"), # Adjust facet label position, alignment, and boldness
      legend.position = "top" # Move the legend to the top of the plot
   )
plot

ggsave(plot, filename = "C:/Users/ives435/OneDrive - PNNL/Desktop/GPCR paper/bpe_b2arTD_BGEcompare.svg",
       scale = 2, width = 7, height = 7, units = "in", dpi = 600)

library(plotly)

# Convert ggplot to plotly
interactive_plot <- ggplotly(plot)

# Display the plot
interactive_plot



