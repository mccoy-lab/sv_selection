library(dplyr)
library(stringr)
library(data.table)
library(pbmcapply)

#########################################################################

### Function for filtering LD results to only include LD between SVs and SNPs.

path_to_files <- "ld_rerun_20200112"

new_cols <- c("SNP_chr", "SNP_pos", "SNP", "SV_chr", "SV_pos", "SV")
get_snp_sv_pairs <- function(pop_id) {
  
  ld_df <- fread(paste0(path_to_files, "/", pop_id, ".ld"))
  
  filtered <- ld_df %>%
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

# get names of population from PLINK within file
populations <- read.table("~/work/syan11/sv_selection/pbs/1KGP_within.txt") %>%
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
  # (this is the copy with the highest r^2 value)
  .[!duplicated(SV), ]

fwrite(all_pops_ld,
       "all_pops_ld.txt",
       sep = "\t")
fwrite(all_pops_max_ld,
       "all_pops_max_ld.txt",
       sep = "\t")


#########################################################################

### Further filtering steps to remove low-quality SVs based on genotyping
### rate, HWE, and allele frequency.

# filter to subset to only high-quality SVs
gt_callrate <- fread("genotyping_callrate.txt")
hwe_pvals <- fread("HWE_pvals.txt")
afs <- fread("eichlerSVs_af_allpops_unfolded.frq.strat.gz") %>%
  group_by(SNP) %>%
  summarize(max_af = max(MAF)) %>%
  setDT()

ld_filtered <- semi_join(all_pops_ld,
                         gt_callrate[ call_rate >= 0.5 ],
                         by = "SV") %>%
  semi_join(., hwe_pvals[ non_hwe_counts <= 13 ],
            by = "SV") %>%
  semi_join(., afs[ max_af > 0 ],
            by = c("SV" = "SNP")) %>%
  as.data.table()


#########################################################################

### some analysis

# how many total SVs are left after filtering?
ld_filtered_max <- ld_filtered %>%
  setorder(., SV, -R2) %>%
  .[!duplicated(SV), ]
nrow(ld_filtered_max) # 63593

# how many SVs have max R2 in any population > 0.8 or > 0.5?
nrow(ld_filtered_max[R2 > 0.8]) # 32423
nrow(ld_filtered_max[R2 > 0.5]) # 42336

# calculate max LD of SV with a SNP in each superpopulation
ld_filtered[pop %in% c("ACB", "ASW", "ESN", "GWD", "LWK", "MSL", "YRI"), superpop := "AFR"]
ld_filtered[pop %in% c("CLM", "MXL", "PEL", "PUR"), superpop := "AMR"]
ld_filtered[pop %in% c("CDX", "CHB", "CHS", "JPT", "KHV"), superpop := "EAS"]
ld_filtered[pop %in% c("CEU", "FIN", "GBR", "IBS", "TSI"), superpop := "EUR"]
ld_filtered[pop %in% c("BEB", "GIH", "ITU", "PJL", "STU"), superpop := "SAS"]
# what percent of SVs have max R2 > 0.8 or > 0.5 in each of the five superpopulations?
lapply(list("AFR", "AMR", "EAS", "EUR", "SAS"),
       function(x) {
         # all SVs in superpop with an LD calculation
         all <- length(unique(ld_filtered[superpop == x]$SV))
         # SVs in superpop with R2 with a SNP > 0.8
         # linked <- length(unique(ld_filtered[superpop == x & R2 > 0.8]$SV))
         linked <- length(unique(ld_filtered[superpop == x & R2 > 0.5]$SV))
         return(paste(x, linked/all))
       }
)


#########################################################################

### Plot the distribution of max LD between SVs and SNPs.

# max LD with a SNP by superpopulation
ld_filtered_max_bypop <- lapply(list("AFR", "AMR", "EAS", "EUR", "SAS"),
                                function(x) ld_filtered[superpop == x] %>%
                                  setorder(., SV, -R2) %>%
                                  .[!duplicated(SV), ]) %>%
  rbindlist()

ld_plot <- ggplot(data = ld_filtered_max_bypop,
                  aes(x = R2, fill = superpop)) +
  geom_histogram(binwidth = 0.05) +
  theme_bw() +
  xlab(expression(paste("Max ",r^2," of SV with a SNP"))) +
  ylab("# of SVs") +
  labs(fill = "Superpopulation") +
  facet_grid(~ superpop, space = "free", scales = "free_x") +
  theme(legend.position = "none",
        panel.grid = element_blank())

ggsave("ld_hist_bypop.pdf",
       ld_plot,
       width = 9,
       height = 2.5)