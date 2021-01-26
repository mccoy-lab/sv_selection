library(data.table)
library(tidyverse)
library(ggplot2)

### Plot a Manhattan-style plot of likelihood ratio statistic (LRS) at the IGH locus,
### including the two adaptive SVs and all SNPs in a 1MB window.
### Color plot background to reflect haplotypes based on the LD heatmap of the region.

############################################################################

### DATA

# path to genotyping callrate results
gt_callrate_path <- "genotyping_callrate.txt"
# path to Hardy-Weinberg p-values results
hwe_pvals_path <- "HWE_pvals.txt"
# path to SV AF results, stratified by population
afs_path <- "eichlerSVs_af_allpops_unfolded.frq.strat.gz"

# path to SV VCF
sv_vcf_path <- "eichlerSVs_1KGP_pgGTs_noseq.vcf.gz"
# path to directory with SV selscan results
sv_selscan_path <- "sv_selscan/"

# path to directory with VCF for SNPs that selscan was run on
snp_vcf_path <- "vcfs/"
# path to directory with selscan output files for SNPs
snp_selscan_path <- "snp_selscan/"

############################################################################


### Read in selscan results for SVs and SNPs

# metadata for filtering SVs
gt_callrate <- fread(gt_callrate_path)
hwe_pvals <- fread(hwe_pvals_path)
afs <- fread(cmd = paste0("zcat ", afs_path))
common_svs <- unique(afs[MAF > 0.05 & MAF < 0.95]$SNP) # common in any 1KGP population
svs <- fread(cmd = paste0("zcat ", sv_vcf_path, " | grep -v '^##' | cut -f 1-3"))

# read in SV selscan results
get_sv_selscan <- function(p) {
  selscan <- fread(paste0(sv_selscan_path, "selscan_50_k8_p", p-1, ".out"))
  # add variant IDs to selscan results
  selscan_svs <- cbind(svs, selscan)
  
  # filter out SVs that violate HWE or have low callrate
  selscan_filt <- semi_join(selscan_svs,
                            gt_callrate[ call_rate >= 0.5 ],
                            by = c("ID" = "SV")) %>%
    semi_join(., hwe_pvals[ non_hwe_counts <= 13 ],
              by = c("ID" = "SV")) %>%
    as.data.table()
  # filter out SVs that are not locally common
  selscan_filt <- selscan_filt[ID %in% common_svs]
  # filter out SVs with extreme global LLE
  selscan_filt <- selscan_filt[`global-lle` > -1000]
  
  setnames(selscan_filt, "lle-ratio", "lle_ratio")
  return(selscan_filt)
}
# get selscan results for ancestry component 2, in which the IGH SVs were under selection
sv_p2 <- get_sv_selscan(2)

# read in SNP selscan results
get_snp_selscan <- function(chr, pop) {
  # read in SNP vcf (for variant IDs) for that chr and ancestry component
  snp_vcf <- fread(cmd = paste0("zcat ", snp_vcf_path,
                                "chr", chr, ".vcf.gz | grep -v '^##' | cut -f 1-3"))
  # read in selscan results for that chr and ancestry component
  selscan <- fread(paste0(snp_selscan_path, "chr", chr,
                          "_selscan_k_p", pop, ".out"))
  
  selscan <- cbind(snp_vcf, selscan)
  setnames(selscan, "lle-ratio", "lle_ratio")
  return(selscan)
}
# get selscan results for ancestry component 2 and chr14, where the IGH SVs are located
snps_chr14_p2 <- get_snp_selscan(14, 2)


############################################################################

### Plot locus-zoom style Manhattan plot for IGH locus

get_locuszoom <- function(sv_ids, pop_number, sv_results_in, snp_results_in, window_size) {
  
  selscan_subset <- sv_results_in[ID %in% sv_ids] %>%
    .[, is_sv := TRUE]
  position <- selscan_subset$POS[1]
  selscan_snp <- snp_results_in %>%
    .[POS > (position - window_size / 2) & POS < (position + window_size / 2)] %>%
    .[, is_sv := FALSE]
  
  selscan_local <- rbind(selscan_subset, selscan_snp) %>%
    .[, locus_title := paste0("Ancestry component ", pop_number)]
  
  return(selscan_local)
}

locuszoom <- get_locuszoom(c("22237_HG02059_ins", "22231_HG02059_del"),
                           2,
                           sv_p2,
                           snp_chr14_p2,
                           1e6)

# sparse points at LRS = 0 for plotting efficiency
locuszoom <- rbind(locuszoom[lle_ratio == 0][sample(1:nrow(locuszoom[lle_ratio == 0]), 1e4)],
                   locuszoom[lle_ratio != 0])

# coordinates for coloring the four haplotypes at the IGH locus,
# based on coordinates for SNPs at edges of the haplotype
haps <- data.frame(start = c(105500865/1000, 105592244/1000, 105711502/1000, 105801600/1000),
                   end = c(105592244/1000, 105711502/1000, 105801600/1000, 105818810/1000),
                   hap = c("A", "B", "C", "D"))

# plot locus zoom plot with haplotype colors
locusplot <- ggplot(data = locuszoom[is_sv == FALSE],
                    aes(x = POS / 1000, y = lle_ratio)) + 
  scale_color_manual(values = "black") +
  # plot haplotype rectangles
  geom_rect(data = haps, aes(NULL, NULL, xmin = start, xmax = end, fill = hap),
            color = "#9E9E9E", size = 0.25, ymin = 0, ymax = Inf, alpha = 0.7) +
  scale_fill_manual(values=c("A" = "#E5A11C",
                             "B" = "#5CB6EA",
                             "C" = "#CC7DAA",
                             "D" = "#2B9F78")) +
  # plot SNPs as black dots
  geom_point(size = 0.75) +
  # plot two SVs as red dots
  geom_point(data = locuszoom[is_sv == TRUE], color = "red", size = 2) +
  ylab("Likelihood ratio statistic") +
  xlab("Position (Kbp)") +
  theme_bw() +
  theme(legend.position = "none", 
        panel.grid = element_blank())

ggsave("locusplot_haplotypes.pdf",
       locusplot,
       width = 3.25, height = 3)