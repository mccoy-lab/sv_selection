library(tidyverse)
library(stringr)
library(data.table)
library(pbmcapply)

setwd("/work-zfs/rmccoy22/syan11/sv_selection/snp_ld")

#########################################################################
### Function that filters LD results to only include LD between SVs and SNPs.
#########################################################################

path_to_files <- "pop_split"
new_cols <- c("SNP_chr", "SNP_pos", "SNP", "SV_chr", "SV_pos", "SV")

get_snp_sv_pairs <- function(pop_id) {
  
  # ld_df <- read.table(paste0(pop_id, "_test.ld"),
  ld_df <- read.table(paste0(path_to_files, "/", pop_id, ".ld"),
                      header = 1,
                      colClasses = c("character", "integer", "character", "character", "integer", "character", "numeric"),
                      stringsAsFactors = FALSE)
  
  filtered <- ld_df %>%
    setDT() %>%
    # filter to only include combinations of a SNP and SV
    .[ ((nchar(SNP_A) > 1) & (nchar(SNP_B) == 1)) |
         ((nchar(SNP_A) == 1) & (nchar(SNP_B) > 1)), ] %>%
    # create new columns for separating SNPs and SVs
    .[, eval(new_cols) := "."] %>%
    # this is possibly the least efficient method known to man
    # assign SV and SNP IDs, positions, etc. to the appropriate columns
    transform(.,
              SNP_pos = ifelse(nchar(SNP_A) == 1,
                               BP_A,
                               BP_B)) %>%
    transform(.,
              SNP_chr = ifelse(nchar(SNP_A) == 1,
                               CHR_A,
                               CHR_B)) %>%
    transform(.,
              SV_pos = ifelse(nchar(SNP_A) > 1,
                              BP_A,
                              BP_B)) %>%
    transform(.,
              SV_chr = ifelse(nchar(SNP_A) > 1,
                              CHR_A,
                              CHR_B)) %>%
    transform(.,
              SV = ifelse(nchar(SNP_A) > 1,
                          SNP_A,
                          SNP_B)) %>%
    # subset columns of interest
    .[, c("SNP_chr", "SNP_pos", "SNP", "SV_chr", "SV_pos", "SV", "R2")]
  
  # sort by SV, then by r^2 value (decreasing)
  unique_svs <- setorder(filtered, SV, -R2) %>%
    # remove rows with duplicate SV names, keeping only the first copy that appears
    # (this is the copy with the highest r^2 value in that population)
    .[!duplicated(SV), ] %>%
    .[, pop := pop_id]
  
  return(unique_svs)
  
}


#########################################################################
### Apply filtering function to all populations. Concatenate the resulting
### dataframes into one, to get an SV's max LD with a SNP in any population.
#########################################################################

# get filepaths of LD files from PLINK
populations <- read.table("~/data/syan11/sv_selection/pbs/1KGP_within.txt") %>%
  .[, c("V3")] %>%
  as.character() %>%
  unique()

# apply filtering function to all populations and concatenate data tables
all_pops_ld <- pbmclapply(populations,
                          get_snp_sv_pairs,
                          mc.cores = getOption("mc.cores", 24L)) %>%
  # concatenate all population data tables into one
  rbindlist()

# get max LD of an SV with a SNP in any of the 26 populations
all_pops_max_ld <- all_pops_ld %>%
  # sort by SV, then by r^2 value (decreasing)
  setorder(., SV, -R2) %>%
  # remove rows with duplicate SV names, keeping only the first copy that appears
  .[!duplicated(SV), ]

#########################################################################
### Plot the distribution of max LD between SVs and SNPs.
#########################################################################

ld_plot <- ggplot(data = all_pops_max_ld,
                  aes(x = R2)) +
  geom_histogram(binwidth = 0.01) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5)) +
  ggtitle("Max LD between SVs and SNPs in all 1KGP populations") +
  xlab("R^2") + ylab("# of SVs")

ggsave("ld_hist.pdf",
       ld_plot,
       width = 9,
       height = 5)