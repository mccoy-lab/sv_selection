library(data.table)
library(tidyverse)
library(ggplot2)
library(rtracklayer)
library(GenomicRanges)

setwd("/scratch/groups/rmccoy22/syan11/sv_selection/ohana")

############################################################################

### Filter selscan results

pbs_results <- fread("pbs_results.txt")
gt_callrate <- fread("genotyping_callrate.txt")
hwe_pvals <- fread("HWE_pvals.txt")
all_pops_ld <- fread("~/work/syan11/sv_selection/snp_ld/all_pops_ld.txt") %>%
  .[, c("SV", "R2", "pop")]
svs <- fread("paragraph_sv_lengths.vcf")

filter_selscan <- function(p) {
  selscan <- fread(paste0("k8/selscan/selscan_50_k8_p",p-1,".out"))
  # selscan_svs <- lazy_dt(cbind(svs, selscan))
  selscan_svs <- cbind(svs, selscan)
  
  # filter out SVs that violate HWE or have low callrate
  selscan_filt <- semi_join(selscan_svs,
                            gt_callrate[ call_rate >= 0.5 ],
                            by = c("ID" = "SV")) %>%
    semi_join(., hwe_pvals[ non_hwe_counts <= 13 ],
              by = c("ID" = "SV")) %>%
    as.data.table()
  
  # filter out SVs where # of steps was > 99
  selscan_filt <- subset(selscan_filt, step < 75)
  
  setorder(selscan_filt, -`lle-ratio`)
  setnames(selscan_filt, "lle-ratio", "lle_ratio")
  # setnames(selscan_filt, "lle-ratio", paste0("lle_ratio_p", p))
  # setnames(selscan_filt, "step", paste0("step_p", p))
  return(selscan_filt)
}

selscan_res <- lapply(1:8, filter_selscan)
# selscan_p2 <- selscan_res[[2]]

############################################################################

### annotations: LD with 1KGP SNPs, genes from gencode

### LD with SNPs
# add superpop ID columns to LD table
afr <- c("YRI", "LWK", "GWD", "MSL", "ESN", "ASW", "ACB")
amr <- c("MXL", "PUR", "CLM", "PEL")
eas <- c("CHB", "JPT", "CHS", "CDX", "KHV")
eur <- c("CEU", "TSI", "FIN", "GBR", "IBS")
sas <- c("GIH", "PJL", "BEB", "STU", "ITU")
all_pops_list <- list(afr, amr, eas, eur, sas) # all populations, grouped by superpop
all_superpops <- c("AFR", "AMR", "EAS", "EUR", "SAS")
all_pops_ld[, superpop := sapply(pop,
                                 get_superpop,
                                 USE.NAMES = FALSE)]

### gene annotations from gencode
# sv metadata
vcf <- fread(cmd = "zcat eichlerSVs_1KGP_pgGTs_noseq.vcf.gz | grep -v '^##' | cut -f 1-3")
vcf <- merge(vcf, svs[, c("ID", "INFO")],
             by = "ID") %>%
  .[, start := POS] %>%
  # calculate end position
  .[, end := as.numeric(NA)] %>%
  .[INFO < 0, end := POS - INFO] %>%
  .[INFO > 0, end := POS]
# convert to granges object
vcf_gr <- makeGRangesFromDataFrame(vcf,
                                   seqnames.field = "#CHROM",
                                   start.field = "start",
                                   end.field = "end",
                                   keep.extra.columns = TRUE)

# gencode gene annotations
gencode <- data.table(readGFF("ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_35/gencode.v35.annotation.gtf.gz"))
genes <- gencode[type == "gene"]
genes_gr <- makeGRangesFromDataFrame(genes, 
                                     seqnames.field = "seqid", 
                                     start.field = "start", 
                                     end.field = "end",
                                     keep.extra.columns = TRUE)

# find overlaps between SVs and gene annotations
olaps <- findOverlaps(vcf_gr, genes_gr)
vcf_matched <- vcf_gr[queryHits(olaps)]
mcols(vcf_matched) <- cbind.data.frame(mcols(vcf_matched),
                                       mcols(genes_gr[subjectHits(olaps)]))
vcf_matched <- vcf_matched %>%
  as.data.table()
sv_gene_overlaps <- group_by(vcf_matched, ID) %>%
  summarize(., genes = list(gene_name)) %>%
  as.data.table()
sv_gene_overlaps[, genes_string := unlist(lapply(sv_gene_overlaps$genes,
                                                 function(x) paste(x, collapse = ", ")))]

ld_genes_annot <- function(dt, superpop_group, i) {
  # get max LD of SV with a SNP in superpop of interest
  ld_subset <- all_pops_ld[superpop == superpop_group,] %>%
    # setorder(., SV, -R2) %>%
    setnames(., "R2", "max_R2_in_superpop") %>%
    # remove rows with duplicate SV names, keeping only the first copy that appears
    # (this is the copy with the highest r^2 value in that population)
    .[!duplicated(SV), ]
  
  # merge selscan results and LD data
  dt <- merge(dt, ld_subset,
              by.x = "ID", by.y = "SV", all.x = TRUE)
  # merge selscan results and gene annotations
  dt <- merge(dt, sv_gene_overlaps,
              by = "ID", all.x = TRUE)
  # add column specifying ancestry component number
  dt$ancestry_component <- paste0("Ancestry component ", i)
  # reorder by highest LLR
  setorder(dt, -lle_ratio)

  return(dt)
}

ancestry_to_superpop <- c("AMR","EAS","EUR","AMR","SAS","EAS","AFR","EUR")
selscan_res2 <- lapply(1:8, function(i) ld_genes_annot(selscan_res[[i]],
                                                       ancestry_to_superpop[i],
                                                       i))
# selscan_p2 <- ld_genes_annot(selscan_p2, "EAS", 2)

# keep only top 10 gene annotations for plotting
less_gene_annot <- function(dt) {
  # only first 10 gene annotations
  keep <- which(!is.na(dt$genes_string))[1:10]
  # set all other gene annotations to NA
  dt[!keep, genes_string := NA]
}
selscan_res2 <- lapply(selscan_res2, function(x) less_gene_annot(x))

############################################################################

### plot of SV length vs. LLR

selscan_all <- rbindlist(selscan_res2)

# no R^2 coloring for SVs with low PBS
selscan_all[lle_ratio < 32, max_R2_in_superpop := NA]

# no gene annotation for SVs with low PBS
# selscan_test[lle_ratio < 80, genes_string := ""]
selscan_all[is.na(genes_string), genes_string := ""]

# pdf(file = "k8/selscan/20201220_p2_R2.pdf", width = 10, height = 5)
pdf(file = "k8/selscan/20201221_p2_R2.pdf", width = 13, height = 10)
ggplot(data = selscan_all,
       # ggplot(data = selscan_plotting,
       aes(x = abs(INFO),
           y = lle_ratio,
           # group = ancestry_component,
           label = genes_string,
           color = max_R2_in_superpop)) +
  geom_point(alpha = 0.3) +
  facet_wrap(~ ancestry_component, scales = "free", nrow = 3) +
  scale_x_log10() +
  theme_classic() +
  xlab("\nlog(SV length (bp))") + ylab("Log likelihood ratio\n") +
  geom_text_repel(size = 3, fontface = "italic") +
  scale_color_viridis_c(option = "plasma", name = "Max R^2 with a known SNP in superpopulation") +
  theme(legend.position = "none") +
  theme(panel.spacing = unit(1, "lines"))
        # strip.background = element_blank(), strip.text = element_text(size = 10))
dev.off()