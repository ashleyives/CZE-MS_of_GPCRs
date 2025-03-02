library(ggplot2)
library(zoo)
library(dplyr)
library(rawrr)
library(pracma)

# File paths
file_paths <- c(
   "C:/Users/ives435/OneDrive - PNNL/Desktop/GPCR paper/massive/20210515_kjb9186_CEMS_mGluR2_30kV_1psi_6.7nL_01.raw",
   "C:/Users/ives435/OneDrive - PNNL/Desktop/GPCR paper/massive/to_add/20210515_kjb9186_CEMS_mGluR2_30kV_1psi_6.7nL_02.raw",
   "C:/Users/ives435/OneDrive - PNNL/Desktop/GPCR paper/massive/to_add/20210515_kjb9186_CEMS_mGluR2_30kV_1psi_6.7nL_03.raw"
)

# Define mass groups as a list
mass_groups <- list(
   groupA = c(1055.27),    # unmodified, 4+
   groupB = c(1075.25),    # monophos, 4+
   groupC = c(1095.24)    # diphos, 4+
)

ppmtol <- 10

# Names of the groups
group_names <- names(mass_groups)

# Initialize an empty data frame to store results
combined_data <- data.frame()

# Process each file
for (file in file_paths) {
   # Initialize an empty list to store chromatogram data for all groups
   C_list <- list()
   
   # Read chromatogram data for each mass group
   for (i in 1:length(mass_groups)) {
      C <- rawrr::readChromatogram(file, mass = mass_groups[[i]], tol = ppmtol, type = "xic", filter = "ms")
      C_list[[i]] <- C
   }
   
   # Process each chromatogram data list
   for (i in 1:length(C_list)) {
      C <- C_list[[i]]
      
      # Extract times and intensities
      times_list <- lapply(C, function(d) d$times)
      intensities_list <- lapply(C, function(d) d$intensities)
      
      # Check if times are identical
      are_identical <- all(sapply(times_list, function(x) identical(x, times_list[[1]])))
      
      # Averaging intensities
      if (are_identical) {
         intensities_matrix <- do.call(cbind, intensities_list)
         y_avg <- rowMeans(intensities_matrix)
         times <- times_list[[1]]
      } else {
         common_times <- sort(unique(unlist(times_list)))
         interpolated_intensities <- sapply(C, function(d) approx(d$times, d$intensities, xout = common_times)$y)
         y_avg <- rowMeans(interpolated_intensities, na.rm = TRUE)
         times <- common_times
      }
      
      # Create data frame and add group and file columns
      data <- data.frame(times = times, y_avg = y_avg, group = group_names[i], file = basename(file))
      
      # Append to combined data
      combined_data <- rbind(combined_data, data)
   }
}

# Print combined data frame
print(combined_data)

# Plot raw averaged values using ggplot2
ggplot(combined_data, aes(x = times, y = y_avg, color = group)) +
   geom_line() +
   labs(title = "Averaged Intensities", y = "Average Intensity", x = "Time") +
   theme_minimal() +
   facet_wrap(~file)

# Gaussian smoothing for all groups
sigma <- 1  # Standard deviation of the Gaussian
window_size <- round(6 * sigma)  # Define the window size, typically 6 * sigma
gaussian_weights <- dnorm(seq(-3, 3, length.out = window_size), mean = 0, sd = sigma)

# Apply rolling window with Gaussian weights for smoothing
combined_data$y_smooth <- NA
for (grp in unique(combined_data$group)) {
   for (fl in unique(combined_data$file)) {
      combined_data$y_smooth[combined_data$group == grp & combined_data$file == fl] <- rollapply(combined_data$y_avg[combined_data$group == grp & combined_data$file == fl], width = window_size, FUN = function(x) sum(x * gaussian_weights), fill = NA, align = "center")
   }
}

# Normalize each group's smoothed intensities to scale to the same maximum
combined_data <- combined_data %>%
   group_by(group, file) %>%
   mutate(y_smooth_scaled = y_smooth / max(y_smooth, na.rm = TRUE)) %>%
   ungroup()

# Plot Gaussian smoothed values with scaling using ggplot2
# Convert 'file' column to a factor
combined_data$file <- as.factor(combined_data$file)

# Check the levels of your 'file' factor
levels(combined_data$file)

# Create a named vector for your labels
facet_labels <- c(`20210515_kjb9186_CEMS_mGluR2_30kV_1psi_6.7nL_01.raw` = "A - Technical Replicate 1", 
                  `20210515_kjb9186_CEMS_mGluR2_30kV_1psi_6.7nL_02.raw` = "B - Technical Replicate 2", 
                  `20210515_kjb9186_CEMS_mGluR2_30kV_1psi_6.7nL_03.raw` = "C - Technical Replicate 3") 

# Create a named vector for legend labels with empty strings to hide them
legend_labels <- c("groupA" = "Unmodified",
                   "groupB" = "Monophosphorylated",
                   "groupC" = "Diphosphorylated")

# Plot
plot <- ggplot(combined_data, aes(x = times, y = y_smooth_scaled, color = group)) +
   geom_line(size = 1) +
   labs(y = "Rel. Abundance", x = "Migration time (min)", color = "mGluR2 Proteoform, 4+") + # Change legend title here
   theme_classic(base_size = 24) +
   xlim(20, 40)+
   coord_cartesian(xlim = c(20, 40))+
   scale_x_continuous(breaks = seq(20, 40, by = 2.5), labels = c(20, "", 25, "", 30, "", 35, "", 40)) +
   facet_wrap(~file, ncol = 1, labeller = as_labeller(facet_labels)) +
   scale_color_manual(values = c("groupA" = "darkblue", "groupB" = "red", "groupC" = "darkgreen"), labels = legend_labels) + # Customize legend labels and colors
   theme(
      text = element_text(family = "Arial"), # Set font to Arial
      strip.text = element_text(face = "bold", family = "Arial", size = rel(1.2)), # Make facet label text bold and bigger
      strip.text.y.left = element_text(angle = 0, hjust = 0, face = "bold", family = "Arial"), # Adjust facet label position, alignment, and boldness
      legend.position = "top" # Move the legend to the top of the plot
   )
plot

ggsave(plot, filename = "C:/Users/ives435/OneDrive - PNNL/Desktop/GPCR paper/xie_grm2_trip.png",
       scale = 2, width = 7, height = 7, units = "in", dpi = 600)

#use plotly to find maxima and calculate migration time/RSD 
library(plotly)

# Convert the ggplot to a plotly plot
plotly_plot <- ggplotly(plot)

# Print the plotly plot
plotly_plot
 
# Values
values <- c(23.69, 23.52, 23.67) #unmodified 
# values <- c(28.59, 28.25, 28.58) #monophos 1 
# values <- c(29.00, 28.71, 28.99) #monophos 2
# values <- c(36.8, 36.4, 37.05) #diphos 1 
# values <- c(37.48, 36.95, 37.70) #diphos 2

# Mean
mean_value <- mean(values)

# Standard Deviation
sd_value <- sd(values)

# RSD %
rsd_percentage <- (sd_value / mean_value) * 100

# Print RSD %
rsd_percentage


