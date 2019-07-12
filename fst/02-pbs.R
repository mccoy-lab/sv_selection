library(data.table)
library(tidyverse)
library(GenomicRanges)
library(annotables)

setwd("/scratch/groups/rmccoy22/syan11/sv_selection/")

afr_ocn <- fread("Africa_Oceania.fst") %>%
  setnames(., c("chr", "pos", "af.afr", "af.ocn", "fst.afr_ocn")) %>%
  mutate(., variant_id = paste(chr, pos, sep = "_")) %>%
  as.data.table()

afr_sas <- fread("Africa_SouthAsia.fst") %>%
  setnames(., c("chr", "pos", "af.afr", "af.sas", "fst.afr_sas")) %>%
  mutate(., variant_id = paste(chr, pos, sep = "_"))  %>%
  as.data.table()

ocn_sas <- fread("Oceania_SouthAsia.fst") %>%
  setnames(., c("chr", "pos", "af.ocn", "af.sas", "fst.ocn_sas")) %>%
  mutate(., variant_id = paste(chr, pos, sep = "_")) %>%
  as.data.table()

results <- merge(merge(afr_ocn, afr_sas, by = "variant_id"), ocn_sas, by = "variant_id") %>%
  mutate(., pbs = ((-log(1 - fst.ocn_sas)) + (-log(1 - fst.afr_ocn)) - (-log(1 - fst.afr_sas))) / 2) %>%
  as.data.table() %>%
  setorder(., -pbs)

genes_gr <- grch38 %>%
  mutate(., strand = "+") %>%
  mutate(., chr = paste0("chr", chr)) %>%
  makeGRangesFromDataFrame(keep.extra.columns = TRUE)

results_gr <- results %>%
  makeGRangesFromDataFrame(keep.extra.columns = TRUE, start.field = "pos", end.field = "pos")

results[, gene_index := nearest(results_gr, genes_gr)]

grch38 <- grch38 %>%
  as.data.table()
grch38[, gene_index := .I]

results <- merge(results, grch38, by = "gene_index", all.x = TRUE) %>%
  setorder(., -pbs)
