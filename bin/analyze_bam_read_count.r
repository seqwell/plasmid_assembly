#!/usr/bin/env Rscript

# Load necessary libraries
library(dplyr)
library(tidyr)
library(purrr)
library(readr)
library(stringr)

# Get arguments: sample_id and bam_readcount file path
args <- commandArgs(trailingOnly = TRUE)
sample_id <- args[1]
bam_file <- args[2]

# Read in the bam-readcount output file
bam_data <- readr::read_table(bam_file, col_names  = FALSE)

# Rename the columns for clarity
colnames(bam_data) <- c("Index", "Position", "RefBase", "TotalReads", "OtherColumns", "BaseCall_A", "BaseCall_C", "BaseCall_G", "BaseCall_T", "BaseCall_N")

# Reshape and process the base calls using pivot_longer instead of gather
df <- bam_data %>%
  select(Position, RefBase, TotalReads, BaseCall_A, BaseCall_C, BaseCall_G, BaseCall_T, BaseCall_N) %>%
  pivot_longer(cols = BaseCall_A:BaseCall_N, names_to = "BaseCallType", values_to = "BaseCallData") %>%
  mutate(counts = str_extract_all(BaseCallData, "[A-Z]:[0-9]+")) %>%
  mutate(counts = str_replace_all(counts, "(A|T|C|G|N):", "")) %>%
  mutate(counts = as.numeric(counts)) 


# Extract base call information
df_clean <- df %>%
  separate(BaseCallData, into = c("Base", "Count"), sep = ":", extra = "drop") %>%
  mutate(Count = as.integer(Count)) %>%
  select(Position, RefBase, TotalReads, BaseCallType, Count)

# Reshape to have separate columns for each base
df_wide <- df_clean %>%
  pivot_wider(names_from = BaseCallType, values_from = Count, values_fill = 0)

# Rename columns for clarity
df_final <- df_wide %>%
  rename(A = `BaseCall_A`,
         C = `BaseCall_C`,
         G = `BaseCall_G`,
         T = `BaseCall_T`,
         N = `BaseCall_N`)


df_final <- df_final %>%
  rowwise() %>%
  mutate(
    Matches = get(RefBase),  # Get the count of matching reads
    Mismatches = TotalReads - Matches  # Compute mismatches
  ) %>%
  ungroup()

# Define output file based on sample_id
output_file <- paste0("per_base_data_", sample_id, ".csv")

# Save the result to a file
readr::write_csv(df_final, output_file)


