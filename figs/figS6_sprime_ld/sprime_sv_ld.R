library(data.table)
library(tidyverse)
library(pbmcapply)
library(ggplot2)
library(rtracklayer)
library(GenomicRanges)
library(liftOver)

### Calculate LD between top hit SVs and SNPs called as introgressed by SPrime,
### to identify adaptive SVs that may also be introgressed.

############################################################################

### DATA

# SPrime introgressed SNP calls: from
# https://md-datasets-cache-zipfiles-prod.s3.eu-west-1.amazonaws.com/y7hyt83vxr-1.zip
sprime_results_path <- "mendeley_data"

# 1000 Genomes sample metadata: from the first sheet of this Excel spreadsheet
# http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/working/20130606_sample_info/20130606_sample_info.xlsx
metadata_path <- "20130606_1KGP_sample_info.csv"

# hg19 to hg38 liftover chain file: from
# http://hgdownload.cse.ucsc.edu/goldenpath/hg19/liftOver/hg19ToHg38.over.chain.gz
chain_path <- "hg19ToHg38.over.chain"

# VCF of 1000 Genomes SNP and SV genotypes
# (generated in this study, in snp_ld section)
snp_sv_concat_path <- "SNP_SV_concat/"

# table of selscan top hit SVs
# (generated in this study)
selscan_results_path <- "selscan_res2.tsv"

# VCF of Paragraph genotypes
# (generated in this study)
vcf_path <- "eichlerSVs_1KGP_pgGTs_noseq.vcf.gz"

############################################################################


# ancestry components assigned to populations
# 1: PUR
# 2: CDX
# 3: FIN
# 4: PEL
# 5: STU
# 6: JPT
# 7: NA (Africa, which should have no introgression)
# 8: TSI
ancestry_comps <- c("PUR", "CDX", "FIN", "PEL", "STU", "JPT", NA, "TSI")

# for a population, collect the sprime output files
get_sprime <- function(pop) {
  # concatenate separate chromosome files
  dt <- rbindlist(lapply(1:22,
                         function(x) fread(paste0(sprime_results_path,
                                                  pop, ".chr", x, ".ND_match"))))
  # add column to mark population
  dt[, population := pop]
  return(dt)
}
# sprime SNPs for each of the 7 ancestry components
sprime <- rbindlist(lapply(as.list(ancestry_comps[!is.na(ancestry_comps)]),
                           function(x) get_sprime(x)))
# add chr prefix to chromosome number
sprime[, CHROM := paste0("chr", CHROM)]
sprime_gr <- makeGRangesFromDataFrame(sprime, 
                                      seqnames.field = "CHROM",
                                      start.field = "POS",
                                      end.field = "POS",
                                      keep.extra.columns = TRUE)
# lift over hg19 coordinates to hg38
sprime_hg38 <- liftOver(sprime_gr,
                        import.chain(chain_path))
sprime_hg38 <- as.data.table(sprime_hg38)[, -c("group", "group_name")]

# make files with list of samples in each 1KGP population of interest (for LD calculation with PLINK)
metadata <- fread(metadata_path)
write_sample_lists <- function(pop) {
  subset <- metadata[Population == pop]
  fwrite(subset[, c("Sample", "Sample")],
         paste0(pop, "_samples.txt"),
         col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")
}
lapply(as.list(ancestry_comps[!is.na(ancestry_comps)]),
       function(x) write_sample_lists(x))

# function to get sprime SNP in highest LD with an SV
get_sv_ld <- function(svID, i, window_size) {
  # name of 1KGP population corresponding to this ancestry component
  i <- as.numeric(i)
  pop <- ancestry_comps[i]
  
  # sprime data for this population
  sprime_data <- sprime_hg38[population == pop]
  
  # SV coordinates
  sv <- top_svs[ID == svID & ancestry_component == paste("Ancestry component", i)]
  chr <- sv[, `#CHROM`]
  pos <- sv[, POS]
  
  # interval to select with tabix
  interval <- paste0(chr, ":", (pos - window_size / 2), "-", (pos + window_size / 2))
  # use tabix to select window around SV from SV/SNP concat VCF
  system(command = paste("tabix",
                         "-h",
                         paste0(snp_sv_concat_path, chr, "_SNP_SV_concat_sort.vcf.gz"),
                         interval,
                         ">",
                         paste0("ld_output/", svID, ".vcf")), intern = TRUE)
  
  # use plink to calculate LD between SV and all other variants in window
  system(command = paste("plink",
                         "--vcf", paste0("ld_output/", svID, ".vcf"),
                         "--keep", paste0(pop, "_samples.txt"),
                         "--r2", 
                         "--ld-snp", svID,
                         "--ld-window 100000",
                         # "--ld-window-r2 0",
                         "--out", paste0("ld_output/", svID)),
         intern = TRUE)
  
  # merge LD output file with sprime SNPs
  ld_out <- fread(paste0("ld_output/", svID, ".ld"))
  ld_out[, CHR_B := paste0("chr", CHR_B)]
  merged <- merge(ld_out, sprime_data,
                  by.x = c("CHR_B", "BP_B"),
                  by.y = c("seqnames", "start"))
  
  # add a dummy row to incorporate information about SV ID, if the merged dt is empty
  dummy_row <- as.list(c(rep(NA, 4), svID, NA, 0, rep(NA, 12)))
  merged <- rbind(merged, dummy_row)
  setorder(merged, -R2)
  
  return(merged[1, ]) # return top LD SNP
}

# wrapper function for `get_sv_ld`
sv_ld_wrapper <- function(row_index, window_size, df) {
  # get SV ID
  svID <- df[row_index, ID]
  print(svID)
  # get ancestry component number
  ancestry_comp_num <- df[row_index, ancestry_component]
  ancestry_comp_num <- gsub("Ancestry component ", "", ancestry_comp_num)
  
  # run function to get max LD of SV with an SPrime introgressed SNP
  output <- get_sv_ld(svID, as.numeric(ancestry_comp_num), window_size)
  return(output)
}

# read in SVs with outlier LRS values
top_svs <- fread(selscan_results_path) %>%
  # throw out SVs significant in ancestry component 7 (Africa, which should have no introgression)
  .[!(ancestry_component == "Ancestry component 7")]

# apply LD wrapper function to every row in `top_svs`
sprime_ld_res <- pbmclapply(1:nrow(top_svs),
                            function(x) sv_ld_wrapper(x, 1e6, top_svs)) %>%
  rbindlist()

# merge with gene annotations
sprime_ld_res <- merge(sprime_ld_res, top_svs[, c("ID", "genes_string")],
                       by.x = c("SNP_A"),
                       by.y = c("ID"),
                       all.x = TRUE)


############################################################################

### filter out SVs that are present in African individuals, who should have no introgression

# make files with list of samples in 1KGP populations of interest,
# for AFR calculation with plink
afr_subset <- metadata[Population %in% c("YRI", "LWK", "GWD", "MSL", "ESN")]
fwrite(afr_subset[, c("Individual ID", "Individual ID")],
       paste0("AFR_nonadmixed_samples.txt"),
       col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")

# calculate AFs with plink
system(command = paste("plink",
                       "--vcf", vcf_path,
                       "--keep AFR_nonadmixed_samples.txt",
                       "--freq gz", 
                       "--keep-allele-order",
                       "--out ld_output/AFR_nonadmixed"),
       intern = TRUE)

# read in plink AFs
afr_afs <- fread("ld_output/AFR_nonadmixed.frq.gz")

# merge with sprime LD data table and only keep SVs absent in Africa (AF < 0.01)
hits_unadmix_afr <- merge(sprime_ld_res, afr_afs[, c("SNP", "MAF")],
                          by.x = c("SNP_A"),
                          by.y = c("SNP"),
                          all.x = TRUE) %>%
  .[MAF < 0.01]

# of unique SVs in high LD (R2 > 0.5) with an SPrime introgressed SNP
length(unique(hits_unadmix_afr[R2 > 0.5]$SNP_A)) # 26


############################################################################

### plotting

# plot distribution of r2 values
p <- ggplot(hits_unadmix_afr, aes(x = as.numeric(R2))) +
  geom_histogram() +
  theme_classic() +
  xlab(expression(paste("Max ", r^2, " with introgressed SPrime SNP"))) +
  ylab("Number of SVs")
ggsave("sprime_r2_hist.pdf", p, width = 5, height = 3)
