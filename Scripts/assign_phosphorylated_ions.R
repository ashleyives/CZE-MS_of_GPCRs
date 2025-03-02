# Load necessary libraries
library(dplyr)
library(tidyverse)
library(Peptides)
library(foreach)
library(doParallel)

#The goal of this script is to take a .csv export from TDValidator make a clean table of which b,y,c, and z ions are Monophosphorylated

#1 - create a list of theoretical fragment ions from a given sequence, this will be used to determine if fragment ions assigned by TDValidator are modified 
########################################################################################################################################################

# Define the sequence of the intatc protein analyte of interest 
#factor xa c terminal cleavage product from b2ar 
sequence <- "TGHGLRRSSKFCLKEHKALKTLGIIMGTFTLCWLPFFIVNIVHVIQDNLIRKEVYILLNWIGYVNSGFNPLIYCRSPDFRIAFQELLCLRRSSLKAYGNGYSSNGNTGEQSGLEVLFQGPYHVEQEKENKLLAEDLPGTEDFVGHQGTVPSDNIDSQGR"

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
#remove internal fragment ions as well 
theoretical_ions_df <- theoretical_ions_df %>%
   filter(!(firstAA == 1 & lastAA == sequence_length))%>%
   filter(firstAA == 1 | lastAA == sequence_length) #remove internals 

# Process the data frame
theoretical_ions_df <- theoretical_ions_df %>%
   # Step 1: Filter out internal entries
   filter(firstAA == 1 | lastAA == sequence_length) %>%
   # Step 2: Duplicate rows with firstAA == 1 and add "b" and "c" types
   mutate(Type = ifelse(firstAA == 1, "B", NA)) %>%
   bind_rows(
      filter(theoretical_ions_df, firstAA == 1) %>%
         mutate(Type = "C")
   ) %>%
   # Step 3: Duplicate rows with lastAA == sequence_length and add "y" and "z" types
   bind_rows(
      filter(theoretical_ions_df, lastAA == sequence_length) %>%
         mutate(Type = "Y")
   ) %>%
   bind_rows(
      filter(theoretical_ions_df, lastAA == sequence_length) %>%
         mutate(Type = "Z")
   ) %>%
   # Remove any additional duplicates if any were created
   distinct() %>%
   filter(!is.na(Type))

# Print the resulting data frame
print(theoretical_ions_df)

# Mutate to create the Type column, adjust monoiso, and create the Name column
theoretical_ions_df <- theoretical_ions_df %>%
   mutate(
      monoiso = case_when(
         Type == "B" ~ monoiso - 18.0105, # water loss 
         Type == "C" ~ monoiso - 1.007, 
         Type == "Z" ~ monoiso - 15.999,
         TRUE ~ monoiso # Default case to retain original monoiso if the type is not B, C, or Z
      ),
      Name = case_when(
         Type %in% c("B", "C") ~ paste(Type, nchar(sequence), sep = ""),
         Type %in% c("Y", "Z") ~ paste(Type, nchar(sequence), sep = ""),
         TRUE ~ NA_character_ # Default case to handle any Type not covered in the conditions
      )
   )

# Print the final data frame
print(theoretical_ions_df)


#2 - import TDValidator export .csv and filter for true positive internal fragment ions 
########################################################################################################################################################

data <- read.csv("C:/Users/ives435/OneDrive - PNNL/Desktop/GPCR paper/20220919_ani0287_cap53_b2ar_pngaseftev30minrt_e2factorxa_et12hcd15cs17_Scan747to772.csv") %>% #monophos 
   left_join(theoretical_ions_df, by = "Name") %>%
   mutate(
      modification_delta = Theoretical.Mass - monoiso,
      modification = if_else(modification_delta >= -1 & modification_delta <= 1, "Unmodified", 
                             if_else(modification_delta >= 79 & modification_delta <= 80, "Monophosphorylated", 
                                     NA_character_))
   ) %>%
   mutate(id = paste(Name, "_", Charge)) 

# Processing the data
ratios <- data %>%
   group_by(id) %>%
   filter(all(c("Monophosphorylated", "Unmodified") %in% modification)) %>% # Keep groups with both modifications
   summarize(ratio = sum(Intensity[modification == "Monophosphorylated"]) / sum(Intensity[modification == "Unmodified"]), .groups = 'drop') %>%
   separate(id, into = c("Name", "Charge"), sep = "_", remove = FALSE) %>%
   separate(Name, into = c("Type", "Position"), sep = "(?<=[A-Za-z])(?=[0-9])", remove = FALSE) %>%
   mutate(Position = as.numeric(Position)) %>%
   mutate(Charge = as.numeric(Charge)) %>%
   mutate(Position = if_else(Type %in% c("Y", "Z"), sequence_length - Position, Position)) 

# Define the positions where vertical lines will be added
vline_positions <- c(8, 9, 92, 93, 107, 111, 151, 156)

ratios %>%
   ggplot(aes(x = Position+268, y = ratio, color = Type, shape=Type)) +
   geom_point(size=3, alpha =0.5) +
   labs(y = "Phosphorylated/Unmodified Fragment Ion Intensity", x = "Amino acid residue") +
   geom_vline(xintercept = vline_positions+268, linetype = "dashed", color = "black") +
   theme_classic(base_family = 24)+
   scale_color_manual(values = c("C" = "black", "Z" = "blue")) # Replace 'Type1' and 'Type2' if your Type variable has specific values

export <- data %>%
   dplyr::select(-Mono..m.z, -id, -firstAA, -lastAA, -sequence, -monoiso, -modification_delta, -prolength, -Type, -Active)

write.csv(export, file="C:/Users/ives435/OneDrive - PNNL/Desktop/GPCR paper/ethcd_factorxa_monophos.csv")


