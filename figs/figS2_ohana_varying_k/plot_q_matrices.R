library(data.table)
library(tidyverse)
library(ggplot2)

### Plot inferred ancestry components for 1000 Genomes samples from Ohana's
### Q matrix, output by `qpas`.

############################################################################

### DATA

# sample IDs from 1000 Genomes, in the same order as in the Q matrix
samples_path <- "samples.txt"

# 1000 Genomes sample metadata: from the first sheet of this Excel spreadsheet
# http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/working/20130606_sample_info/20130606_sample_info.xlsx
metadata_path <- "20130606_1KGP_sample_info.csv"

# Q matrices, output from running `qpas` with different k values

############################################################################


# sample order in the Q matrix
sid <- fread(samples_path, header = FALSE) %>%
  setnames(., "sample_id")

# 1000 Genomes metadata with population assignments for samples
pid <- fread(metadata_path) %>%
  setnames(., "Sample", "sample_id")

plot_q_matrix <- function(matrix_path, k) {
  # read Q matrix from Ohana and reformat
  qmat <- fread(matrix_path, skip = 1) %>%
    .[, sample_id := sid$sample_id] %>%
    pivot_longer(
      cols = starts_with("V"),
      names_to = "ancestry_component",
      names_prefix = "V",
      values_to = "proportion",
      values_drop_na = TRUE
    )
  
  # merge to with 1000 Genomes metadata get sample and population labels
  qmat <- merge(qmat, pid, by = "sample_id") %>%
    setDT()
  
  # reorder based on superpopulation
  qmat$Population <- factor(qmat$Population,
                            # levels = pop_vector)
                            levels = c("ACB","ASW","ESN","GWD","LWK","MSL","YRI", "CDX","CHB","CHS","JPT","KHV",
                                       "CEU","FIN","GBR","IBS","TSI", "BEB","GIH","ITU","PJL","STU",
                                       "CLM","PEL","PUR","MXL"))
  
  # plot
  p <- ggplot(data = qmat, aes(x = sample_id, y = proportion, fill = ancestry_component)) +
    geom_bar(position = "stack", stat = "identity") +
    # separate subplots for each population
    facet_wrap(~ Population, scales = "free", nrow = 1) +
    # ggtitle(paste0("k = ", k)) +
    theme(axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.title.x = element_blank(), axis.title.y = element_blank(),
          strip.text.x = element_text(size = 12, face = "bold"))
  
  ggsave(paste0("qmat_k", k, ".pdf"), p, width = 28, height = 1.5)
}

# plot Q matrix for the different k values we attempted
lapply(as.list(c(5, 8, 12, 16)),
       function(x) plot_q_matrix(paste0("k", x, "/chr21_pruned_50_Q.matrix"),
                                 x))