library(rstudioapi)
library(dplyr)

# Getting the path of your current open file
current_path = rstudioapi::getActiveDocumentContext()$path 
setwd(dirname(current_path ))
getwd()

hmmsearch_tblout <- read.csv("../../results/hmmsearch/hmmsearch_tblout.txt", skip = 0, sep = "-", header = FALSE, strip.white = TRUE)

hmmsearch_tblout <- hmmsearch_tblout %>% filter(!grepl('#', V1))
hmmsearch_tblout <- rename(hmmsearch_tblout[1], pephash = V1)



prodigal_all <- read.csv("../../results/prodigal/prodigal_out.all.nr_expanded.csv")

hmm_hits <- left_join(hmmsearch_tblout, prodigal_all, by = 'pephash')

write.csv(hmm_hits, "../../results/hmmsearch/hmm_hits.csv")
