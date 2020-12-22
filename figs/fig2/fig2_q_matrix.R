#### make ancestry component plots, separated by superpop
#### and ordered by admixture proportion

library(tidyverse)
library(data.table)
library(ggplot2)

setwd("/Users/syan/Documents/mccoy-lab/sv_selection/ohana")

# sample names
sid <- fread("samples.txt", header = FALSE) %>%
  setnames(., "sample_id")
# samples with population ID
pid <- fread("20130606_g1k.ped") %>%
  setnames(., "Individual ID", "sample_id")

qmat <- fread("k8/chr21_pruned_50_Q.matrix", skip = 1) %>%
  .[, sample_id := sid$sample_id] %>%
  pivot_longer(
    cols = starts_with("V"),
    names_to = "ancestry_component",
    names_prefix = "V",
    values_to = "proportion",
    values_drop_na = TRUE
  )

# add in which population a sample belongs to
qmat <- merge(qmat, pid, by = "sample_id") %>%
  as.data.table()

# superpop: string with superpopulation ID
# pop_names: vector of populations in that superpop
# comp_num: ancestry component # to stratify by
reorder_qmat <- function(superpop, pop_names, comp_num) {
  # select only populations of interest
  subset <- qmat[Population %in% pop_names]
  # group samples by population for plotting
  subset$Population <- factor(subset$Population,
                              levels = c(pop_names))
  # reorder samples by ancestry component of interest
  order <- subset[ancestry_component == comp_num] %>%
    setorder(., -proportion) %>%
    .$sample_id
  # apply reordering to dataframe for plotting
  subset$sample_id <- factor(subset$sample_id, levels = order)
  return(subset)
}

afr <- reorder_qmat("afr", c("ACB","ASW","ESN","GWD","LWK","MSL","YRI"), 7)
eas <- reorder_qmat("eas", c("CHB","CHS","CDX","KHV","JPT"), 2)
sas <- reorder_qmat("sas", c("BEB","GIH","ITU","PJL","STU"), 5)
eur <- reorder_qmat("eur", c("CEU","FIN","GBR","IBS","TSI"), 8)
amr <- reorder_qmat("amr", c("CLM","MXL","PEL","PUR"), 8)
qmat_reordered <- rbind(afr, eas, sas, eur, amr)

p <- ggplot(data = qmat_reordered,
            aes(x = sample_id, y = proportion, fill = ancestry_component)) +
  geom_bar(position = "stack", stat = "identity") +
  facet_wrap(~ Population, scales = "free", nrow = 1) +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        strip.text = element_text(size = 15, face = "bold"),
        legend.position = "none")

ggsave(paste0("~/Documents/mccoy-lab/presentations/sv_selection_paper/fig2/20201217_qmat_50_reordered.pdf"),
       p,
       width = 40, height = 2.1)
