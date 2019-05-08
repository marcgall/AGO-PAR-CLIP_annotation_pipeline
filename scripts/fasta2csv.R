#!/usr/bin/env Rscript
#
# Invoke with Rscript --vanilla fasta2csv.R miR.fa miR.csv
# arg1: input file.
# arg2: output file
args = commandArgs(trailingOnly=TRUE)

library(Biostrings)
library(dplyr)
library(readr)

if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)}

fa <- 
  readBStringSet(args[1]) %>%
  as.character() %>% 
  tibble::enframe(name = "id", value = "seq") %>%
  filter(grepl("Homo sapiens", id)) %>%
  mutate(id = regmatches(x = id, #trimming id
                         regexpr(pattern = "^\\S+", 
                                 text = id))) 

if (length(args) == 1) {
  output <- 
    strsplit(args[1], split = ".", fixed = TRUE) %>% 
    .[[1]] %>% .[1] %>% paste0(., ".csv")
} else {
  output <- args[2]}

write_csv(fa, path = output, col_names = FALSE)