library(ggplot2)
library(zoo)
library(dplyr)
library(rawrr)

# File paths
file_paths <- c(
   "C:/Users/ives435/OneDrive - PNNL/Desktop/GPCR paper/massive/20220920_ani0287_cap53_b2ar_pngaseftev30minrt_e2factorxa_et12hcd15_1066p24.raw",
   "C:/Users/ives435/OneDrive - PNNL/Desktop/GPCR paper/massive/20220920_ani0287_cap53_b2ar_pngaseftev30minrt_e2factorxa_et12hcd25_906p41.raw",
   "C:/Users/ives435/OneDrive - PNNL/Desktop/GPCR paper/massive/20220920_ani0287_cap53_b2ar_pngaseftev30minrt_e2factorxa_et12hcd25_954p01.raw"
)

# Initialize an empty data frame to store results
combined_data <- data.frame()

# Process each file
for (file in file_paths) {
   # Read chromatogram data for each file
   C <- rawrr::readChromatogram(file, type = "tic")
   
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
facet_labels <- c(`20220920_ani0287_cap53_b2ar_pngaseftev30minrt_e2factorxa_et12hcd15_1066p24.raw` = "A - Technical Replicate 1", 
                  `20220920_ani0287_cap53_b2ar_pngaseftev30minrt_e2factorxa_et12hcd25_906p41.raw` = "B - Technical Replicate 2", 
                  `20220920_ani0287_cap53_b2ar_pngaseftev30minrt_e2factorxa_et12hcd25_954p01.raw` = "C - Technical Replicate 3") 

# Plot
plot <- ggplot(combined_data, aes(x = times, y = intensities)) +
   geom_line(size = 1) +
   geom_vline(xintercept = c(26, 32), linetype = "dashed", color = "blue", size=2) + # Adding vertical lines
   geom_vline(xintercept = c(10, 25), linetype = "dashed", color = "black", size=2) + # Adding vertical lines
   labs(y = "Rel. Abundance", x = "Migration time (min)") + # Change legend title here
   theme_classic(base_size = 24) +
   ylim(0, 3e9)+
   facet_wrap(~file, ncol = 1, labeller = as_labeller(facet_labels)) +
   theme(
      text = element_text(family = "Arial"), # Set font to Arial
      strip.text = element_text(face = "bold", family = "Arial", size = rel(1.2)), # Make facet label text bold and bigger
      strip.text.y.left = element_text(angle = 0, hjust = 0, face = "bold", family = "Arial"), # Adjust facet label position, alignment, and boldness
      legend.position = "top" # Move the legend to the top of the plot
   )
plot

ggsave(plot, filename = "C:/Users/ives435/OneDrive - PNNL/Desktop/GPCR paper/tic_b2arMD_trip.png",
       scale = 2, width = 7, height = 7, units = "in", dpi = 600)






