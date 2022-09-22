
# Specify packages
my_packages <- c("dplyr")
# Extract not installed packages
not_installed <- my_packages[!(my_packages %in% installed.packages()[ , "Package"])]
# Install not installed packages
if (length(not_installed)) install.packages(not_installed)


#library(rstudioapi)
library(dplyr)

# Getting the path of your current open file
# current_path = rstudioapi::getActiveDocumentContext()$path 
# setwd(dirname(current_path ))
# getwd()

hmmsearch_tblout <- read.csv(snakemake@input[['tblout']], skip = 0, sep = "-", header = FALSE, strip.white = TRUE)

hmmsearch_tblout <- hmmsearch_tblout %>% filter(!grepl('#', V1))
hmmsearch_tblout <- rename(hmmsearch_tblout[1], pephash = V1)

prodigal_all <- read.csv(snakemake@input[['prodigal_all']])

hmm_hits <- left_join(hmmsearch_tblout, prodigal_all, by = 'pephash')

write.csv(hmm_hits, snakemake@output[['out']])
