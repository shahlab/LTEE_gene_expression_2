library(tidyverse)
library(broom)
library(parallel)

# you may need to set working directory
#setwd("code/analysis/")

# load the deseq results for rnaseq
ddf <- read_csv("../../data_frames/table_s2_fold_changes.csv") %>%
  filter(seqtype == "rna")

# make a separate df for each line, removing unnecessary columns
line.dfs <- ddf %>%
  filter(seqtype == "rna") %>%
  select(line, target_id, log2FoldChange, padj) %>%
  split(.$line)

# how many genes were significant in each line, in each direction?
obs.sig.counts <- ddf %>%
  filter(padj <= .01) %>%
  mutate(direc = ifelse(log2FoldChange > 0, "up", "down")) %>%
  group_by(target_id, direc) %>%
  tally() %>%
  ungroup() %>%
  group_by(direc, n) %>%
  tally() %>%
  ungroup() %>%
  add_case(direc = "up", n = 10, nn = 0) # there's no up n=10 case, add it.

# run the randomization a million times
obs.over.exp <- mclapply(1:1e6, function(x) {
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
  
  # create the expected distribution by counting the number of times
  # a gene is sig across the lines for up/downregs
  expected <- shuffled.names %>%
    filter(padj <= .01) %>%
    mutate(direc = ifelse(log2FoldChange > 0, "up", "down")) %>%
    group_by(target_id, direc) %>%
    tally() %>%
    ungroup() %>%
    group_by(direc, n) %>%
    tally(name = "nn") %>%
    ungroup()
  
  # then compare it to the observed dist by doing a join
  full_join(
    obs.sig.counts,
    expected,
    by = c("direc", "n"),
    suffix = c("_obs", "_exp")
  ) %>%
    mutate(ratio = nn_obs / nn_exp) # find the ratio of obs / exp
}, mc.cores = 32) %>%
  bind_rows(.id = "iteration")

count.cis <- obs.over.exp %>%
  # if expected value isn't present because it never occurred, make it 0
  mutate(nn_exp = ifelse(is.na(nn_exp), 0, nn_exp)) %>% 
  # for each direction and number of lines a gene was sig in (n)
  group_by(direc, n) %>%
  summarise(
    mean_count = mean(nn_exp, na.rm = TRUE), # the mean expected value
    ci = 2.576 * (sd(nn_exp, na.rm = TRUE) / sqrt(n())), # the CI
    ci.lower = mean_count - ci,
    ci.upper = mean_count + ci
  ) %>%
  ungroup() %>%
  mutate(
    # mean prob is mean number of expected genes / total genes
    mean_prob = mean_count / length(unique(ddf$target_id)) 
  )

write_csv(count.cis, "/data/john/projects/3ltee/data_frames/count_cis.csv")
