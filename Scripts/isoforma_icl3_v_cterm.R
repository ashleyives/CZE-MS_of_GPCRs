library(isoforma)
library(pspecterlib)
library(dplyr)
library(Rdisop)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(isopat)

file <- "C:/Users/ives435/OneDrive - PNNL/Desktop/GPCR paper/massive/20220919_ani0287_cap53_b2ar_pngaseftev30minrt_e2factorxa_et12hcd15cs17.raw" #monophos 
Scan_Numbers <- c(seq(742, 772, by=2))

xml_data <- get_scan_metadata(MSPath = file)

Sequence <- "TGHGLRRSSKFCLKEHKALKTLGIIMGTFTLCWLPFFIVNIVHVIQDNLIRKEVYILLNWIGYVNSGFNPLIYCRSPDFRIAFQELLCLRRSSLKAYGNGYSSNGNTGEQSGLEVLFQGPYHVEQEKENKLLAEDLPGTEDFVGHQGTVPSDNIDSQGR"

positions <- gregexpr("[STY]", Sequence)[[1]]

# Convert positions to a space-separated string
positions_string <- paste(positions, collapse = " ")

# Print the result
print(positions_string)

Summed_Peaks <- sum_ms2_spectra(ScanMetadata = xml_data, ScanNumbers = Scan_Numbers)

head(Summed_Peaks)

Modifications <- "Phospho,(9, 92)[1]" #divides into icl3 versus c-terminus 
Modified_Sequences <- pspecterlib::multiple_modifications(Sequence, Modifications, ReturnUnmodified = TRUE)

IsoForma_Example <- isoforma_pipeline(
   Sequences = Modified_Sequences,
   SummedSpectra = Summed_Peaks,
   PrecursorCharge = 17, # just a maximum charge threshold so that it doesn't look for fragment ions bigger than this #
   ActivationMethod = "ETD", #if not specificed as HCD or ETD will search all ion types 
   IonGroup = "c",
   IsotopeAlgorithm = "isopat", # Rdisop is preferred, is faster, and is more accurate, but it tends to crash on Windows
   Message = TRUE
)

IsoForma_Example[[1]]
IsoForma_Example[[2]]
IsoForma_Example[[3]]
IsoForma_Example[[4]]
IsoForma_Example[[5]]