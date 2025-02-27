# Load necessary libraries
library(dplyr)
library(tidyverse)
library(Peptides)
library(foreach)
library(doParallel)

#The goal of this script is to take a .csv export from TDValidator (w/ tentative internal fragment ions assigned) and eliminate all false positive internal fragments ions 
#To sanity check if this is working, put a mock/scrambled sequence into TDValidator and ensure that the filtering steps shown below remove all false positive internal fragment ions

#1 - create a list of theoretical fragment ions from a given sequence, this will be used to determine if fragment ions assigned by TDValidator are modified 
########################################################################################################################################################

# Define the sequence of the intatc protein analyte of interest 
sequence <- "DYKDDDDAMGQPGNGSAFLLAPNRSHAPDHDVENLYFQGTQQRDEVWVVGMGIVMSLIVLAIVFGNVLVITAIAKFERLQTVTNYFITSLACADLVMGLAVVPFGAAHILTKTWTFGNFWCEFWTSIDVLCVTASIETLCVIAVDRYFAITSPFKYQSLLTKNKARVIILMVWIVSGLTSFLPIQMHWYRATHQEAINCYAEETCCDFFTNQAYAIASSIVSFYVPLVIMVFVYSRVFQEAKRQLQKIDKSEGRFHVQNLSQVEQDGRTGHGLRRSSKFCLKEHKALKTLGIIMGTFTLCWLPFFIVNIVHVIQDNLIRKEVYILLNWIGYVNSGFNPLIYCRSPDFRIAFQELLCLRRSSLKAYGNGYSSNGNTGEQSGLEVLFQGPYHVEQEKENKLLAEDLPGTEDFVGHQGTVPSDNIDSQGRNASTNDSLLETSQVAPA"

# Calculate the total length of the sequence
sequence_length <- nchar(sequence)

# Create the prolength column value to be used in each row
prolength_value <- sequence_length

# Initialize a data frame to store truncations
theoretical_ions_df <- data.frame(firstAA = numeric(), lastAA = numeric(), sequence = character(), monoiso = numeric(), prolength = numeric(), stringsAsFactors = FALSE)

# Register parallel backend
num_cores <- detectCores() - 1  # Use one less than the total number of cores
cl <- makeCluster(num_cores)
registerDoParallel(cl)

# Parallelized loop over each starting position in the sequence
results <- foreach(start_pos = 1:sequence_length, .combine = rbind, .packages = c("Peptides", "dplyr")) %dopar% {
   # Initialize a temporary data frame to store truncations for each start_pos
   temp_df <- data.frame(firstAA = numeric(), lastAA = numeric(), sequence = character(), monoiso = numeric(), prolength = numeric(), stringsAsFactors = FALSE)
   
   # Loop through the sequence to generate truncations from start_pos to the end
   for (end_pos in start_pos:sequence_length) {
      truncation <- substr(sequence, start_pos, end_pos)
      monoiso_weight <- Peptides::mw(truncation, monoisotopic = TRUE, avgScale = "expasy")
      temp_df <- rbind(temp_df, data.frame(firstAA = start_pos, lastAA = end_pos, sequence = truncation, monoiso = monoiso_weight, prolength = prolength_value, stringsAsFactors = FALSE))
   }
   
   return(temp_df)
}

# Stop the parallel backend
stopCluster(cl)

# Convert the results to a data frame
theoretical_ions_df <- as.data.frame(results)

# Filter out the sequences where firstAA == 1 and lastAA == sequence_length
theoretical_ions_df <- theoretical_ions_df %>%
   filter(!(firstAA == 1 & lastAA == sequence_length))

# Mutate to create the Type column, adjust monoiso, and create the Name column
theoretical_ions_df <- theoretical_ions_df %>%
   mutate(
      Type = case_when(
         firstAA == 1 ~ "B",
         lastAA == sequence_length ~ "Y",
         TRUE ~ "I"
      ),
      monoiso = case_when(
         Type == "B" ~ monoiso - 18.0105, #water loss 
         Type == "I" ~ monoiso - 18.0105, #water loss
         TRUE ~ monoiso
      ),
      Name = case_when(
         Type == "B" ~ paste(Type, nchar(sequence), sep = ""),
         Type == "Y" ~ paste(Type, nchar(sequence), sep = ""),
         Type == "I" ~ paste(Type, firstAA-1, "-", lastAA, sep = "") #shift firstAA for internal ion naming 
      )
   )

# Print the final data frame
print(theoretical_ions_df)

#2 - import TDValidator export .csv and filter for true positive internal fragment ions 
########################################################################################################################################################

data0 <- read.csv("C:/Users/ives435/OneDrive - PNNL/Desktop/GPCR paper/20220915_cap53_b2arpngasef20220818_HCD20_monophos_scan942to1348-qb.csv") %>%
   mutate(Mono.m.z = Mono..m.z,
          roundmz = round(Mono.m.z, 0)) %>%
   mutate(id = paste(roundmz, "_", Charge)) %>%
   dplyr::select(-Active)

#all data is shifted by 2-3 ppm, will re-calibrate ppms below  
hist(data0$ppm.Error, 100)

data <- data0 %>%
   mutate(adj.ppm.Error = ppm.Error-median(ppm.Error))%>%
   separate_rows(Name, sep = ", ") #if 'Name' lists two ions splits into two rows 

#find fragment ions that have same mz, can also do MZ_charge  
#if there is a duplicate, preferentially pick "B" or "Y" over internals "I"
#return only internal fragment ions with |adjusted ppm error| <2 ppm 
filtered_data <- data %>%
   group_by(roundmz) %>%
   filter(if(any(grepl("B|Y", Name))) grepl("B|Y", Name) else grepl("I", Name)) %>% 
   slice_max(order_by = Score, n = 1) %>% #if ions share the name mz, return which ever row has the best isotopic fit score
   filter(ifelse(grepl("B|Y", Name), adj.ppm.Error >= -10 & adj.ppm.Error <= 10, adj.ppm.Error >= -5 & adj.ppm.Error <= 5)) %>%
   filter(Mono..m.z < 1800) %>% #above 1800 noise is very high
   left_join(theoretical_ions_df, by = "Name") %>%
   mutate(modification_delta = Theoretical.Mass-monoiso)

export <- filtered_data %>%
   dplyr::select(-Mono..m.z, -roundmz, -id, -firstAA, -lastAA, -sequence, -monoiso, -prolength, -Type)

write.csv(export, file="C:/Users/ives435/OneDrive - PNNL/Desktop/GPCR paper/filtered_frags_HCD20_monophos.csv")

#sanity if there are ions with duplicate m/z and charge, can happen if two internal fragment ions have same isotopic fit score
duplicated_ids <- filtered_data %>%
   group_by(id) %>%
   filter(n() > 1) %>%
   distinct(id)






