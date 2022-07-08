library(tidyverse)
library(broom)
library(parallel)

# you may need to set working directory
setwd("code/analysis/")

# load the deseq results for rnaseq
ddf <- read_csv("../../data_frames/table_s2_fold_changes.csv") %>%
  filter(seqtype == "rna")

# make a separate df for each line, removing unnecessary columns
line.dfs <- ddf %>%
  filter(seqtype == "rna") %>%
  select(line, target_id, log2FoldChange, padj) %>%
  split(.$line)

# these genes are to be removed, they are genes that are counting towards up and
# down, i.e. those that are up in 4 and down in 4, etc
double.genes <- ddf %>%
  filter(padj <= .01) %>%
  mutate(direc = ifelse(log2FoldChange > 0, "up", "down")) %>%
  group_by(target_id, direc) %>%
  tally() %>%
  ungroup() %>%
  group_by(target_id) %>%
  tally() %>%
  arrange(desc(n)) %>%
  filter(n > 1) %>%
  pull(target_id)

# how many genes were significant in each line, in each direction?
obs.sig.counts <- ddf %>%
  filter(padj <= .01 & !(target_id %in% double.genes)) %>%
  mutate(direc = ifelse(log2FoldChange > 0, "up", "down")) %>%
  group_by(target_id, direc) %>%
  tally() %>%
  ungroup() %>%
  group_by(direc, n) %>%
  tally() %>%
  ungroup() %>%
  add_case(direc = "up", n = 10, nn = 0) # there's no up n=10 case, add it.

# run the randomization a million times
obs.over.exp <- mclapply(1:1e4, function(x) {
  # shuffle the names about the values for each line, combining to 1 df later.
  # shuffling the names about the l2fc and padj is easier because it keeps the
  # l2fc and padj paired.
  shuffled.names <- line.dfs %>%
    map(function(y) {
      # copy the df
      tmp <- y
      
      # shuffle the names about the l2fcs/padj for each df
      tmp$target_id <-
        sample(tmp$target_id, size = nrow(tmp), replace = FALSE)
      
      # return the df with shuffled names
      return(tmp)
    }) %>%
    bind_rows() # combine to a single df
  
  # that randomization produces a new set of double genes, find them again and
  # exclude them
  dg <- shuffled.names %>%
    filter(padj <= .01) %>%
    mutate(direc = ifelse(log2FoldChange > 0, "up", "down")) %>%
    group_by(target_id, direc) %>%
    tally() %>%
    ungroup() %>%
    group_by(target_id) %>%
    tally() %>%
    arrange(desc(n)) %>%
    filter(n > 1) %>%
    pull(target_id)
  
  data.frame(ndub=length(dg))
}, mc.cores = 32) %>%
  bind_rows(.id = "iteration")

write_csv(obs.over.exp, "/data/john/projects/3ltee/data_frames/dg_per_sim.cs")
