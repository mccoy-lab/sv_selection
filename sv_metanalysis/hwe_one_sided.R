library(HardyWeinberg)
library(data.table)
library(dplyr)
library(pbmcapply)

### Calculates one-sided Hardy-Weinberg p-values in each 1000 Genomes population
### and counts the number of populations in which an SV violates one-sided HWE
### (excess of heterozygotes). Y chromosome variants are ignored because they're haploid.

############################################################################


get_hwe_pvals <- function(population_ID) {
  
  autosomal_file <- paste0("gt_counts/",
                           population_ID,
                           "_gtcounts_autosome.txt")
  xchrom_file <- paste0("gt_counts/",
                        population_ID,
                        "_gtcounts_xchr.txt")
  
  # Do HWE calculations for autosomal variants
  
  # data table columns should be ordered: homozygous, het, homozygous
  # order of homozygotes for the minor vs. major allele does not matter
  a_data <- fread(autosomal_file,
                  header = FALSE,
                  col.names = c("SV","AA","AB","BB","missing"))
  
  # convert data into a matrix for HWExactMat
  a_data_mat <- data.matrix(a_data[,c("AA","AB","BB")])
  
  # calculate one-sided HWE p-values,
  # where only an excess of heterozygotes is evidence against HWE
  a_data$pvals <- HWExactMat(a_data_mat,
                             alternative = "greater",
                             x.linked = FALSE,
                             verbose = FALSE)$pvalvec
  
  
  # Do HWE calculations for variants on the X chromosome
  
  # X-chromosome input dataframe has genotype counts in order A, B, AA, AB, BB
  # A and B are the male counts, the rest are the female counts
  x_data <- fread(xchrom_file,
                  header = FALSE,
                  col.names = c("SV","A","B","AA","AB","BB","missing"))

  # convert data into a matrix for HWExactMat
  x_data_mat <- data.matrix(x_data[,c("A","B","AA","AB","BB")])

  # calculate one-sided HWE p-values,
  # where only an excess of heterozygotes is evidence against HWE
  x_data$pvals <- HWExactMat(x_data_mat,
                             alternative = "greater",
                             x.linked = TRUE,
                             verbose = FALSE)$pvalvec


  # combine autosome and X chromosome variants
  data <- rbind(a_data[,c("SV","pvals")],
                x_data[,c("SV","pvals")])
  return(data$pvals)
  
}


# get list of populations
file_list <- list.files(path = "gt_counts",
                        pattern = ".txt")
populations <- sapply(strsplit(file_list, "_"), `[`, 1) %>%
  unique()

# get vector of SV names from one of the gt counts files
autosomal_SVs <- fread("gt_counts/ACB_gtcounts_autosome.txt",
                       header = FALSE)[,1]
xchr_SVs <- fread("gt_counts/ACB_gtcounts_xchr.txt",
                  header = FALSE)[,1]

# get HWE p-values for all 1000 Genomes population
pvals <- pbmclapply(populations,
                    FUN = get_hwe_pvals,
                    mc.cores = getOption("mc.cores", 24L)) %>%
  setDT() %>%
  setnames(populations) %>%
  .[, SV := c(autosomal_SVs, xchr_SVs)]

# count # of populations where SV's HW p-value is < 0.0001
pvals_count <- mutate_all(pvals[,-27],
                          function(x) {ifelse((x < 0.0001),
                                              1,
                                              0)}) %>%
  mutate(non_hwe_counts = rowSums(.)) %>%
  setDT() %>%
  .[, SV := c(autosomal_SVs, xchr_SVs)]

# add column to pvals table with counts of populations in which SV violates HWE
pvals %>%
  .[, non_hwe_counts := pvals_count[, non_hwe_counts]]

fwrite(pvals, "HWE_pvals.txt",
       sep = "\t")