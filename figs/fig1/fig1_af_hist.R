library(data.table)
library(tidyverse)
library(ggplot2)

#########################################################################
### Read in AFs from PLINK and plot histograms by continental population.
#########################################################################

setwd("/work-zfs/rmccoy22/syan11/sv_selection/af_plink/")

gt_callrate_path <- "../paragraph/post_processing/genotyping_callrate.txt"
hwe_path <- "../paragraph/post_processing/HWE_pvals.txt"
af_path <- "eichlerSVs_af_allpops_unfolded.frq.strat.gz"

afs <- fread(af_path)

afs <- semi_join(afs,
                 gt_callrate[ call_rate >= 0.5 ],
                 by = c("SNP" = "SV")) %>%
  semi_join(., hwe_pvals[ non_hwe_counts <= 13 ],
            by = c("SNP" = "SV")) %>%
  as.data.table()

# group SVs by continental population
afr <- c("YRI", "LWK", "GWD", "MSL", "ESN", "ASW", "ACB")
amr <- c("MXL", "PUR", "CLM", "PEL")
eas <- c("CHB", "JPT", "CHS", "CDX", "KHV")
eur <- c("CEU", "TSI", "FIN", "GBR", "IBS")
sas <- c("GIH", "PJL", "BEB", "STU", "ITU")

all_pops_list <- list(afr, amr, eas, eur, sas) # all populations, grouped by superpop
all_superpops <- c("AFR", "AMR", "EAS", "EUR", "SAS")
get_superpop <- function(pop) {
  idx_vector <- sapply(all_pops_list, function(y) pop %in% y)
  return(all_superpops[idx_vector])
}
afs[, superpop := sapply(CLST,
                         get_superpop,
                         USE.NAMES = FALSE)]

afs_eas <- afs[superpop == "EAS",]
afs_afr <- afs[superpop == "AFR",]
afs_amr <- afs[superpop == "AMR",]
afs_eur <- afs[superpop == "EUR",]
afs_sas <- afs[superpop == "SAS",]

# plot distribution of AFs per continental population
af_plot <- ggplot(data = afs_eas,
                  aes(x = MAF,
                      color = superpop)) +
  geom_density(alpha = 0.1,
               # divide to normalize the SV count and make it look like the histogram bar height
               aes(y = (..count..) / 100 )) +
  geom_density(alpha = 0.1,
               data = afs_afr,
               aes(x = MAF,
                   color = superpop,
                   y = (..count..) / 120 )) +
  geom_density(alpha = 0.1,
               data = afs_amr,
               aes(x = MAF,
                   color = superpop,
                   y = (..count..) / 80 )) +
  geom_density(alpha = 0.1,
               data = afs_eur,
               aes(x = MAF,
                   color = superpop,
                   y = (..count..) / 100 )) +
  geom_density(alpha = 0.1,
               data = afs_sas,
               aes(x = MAF,
                   color = superpop,
                   y = (..count..) / 100 )) +
  # histogram with the same bins
  # geom_histogram(alpha = .5, data = afs_sas, aes(x = MAF, color = superpop, y=..count../5)) +
  xlim(0,1) +
  # scale_y_continuous(trans='log2') + # this looks TERRIBLE
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5)) +
  ggtitle("AF for all SVs genotyped in 1KGP") +
  xlab("Allele frequency") + ylab("Number of SVs")

ggsave("afs_filtered_count_0927.pdf",
       af_plot,
       width = 9,
       height = 5)