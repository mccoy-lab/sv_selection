library(data.table)
library(tidyverse)
library(ggplot2)
library(rtracklayer)
library(GenomicRanges)
library(ggrastr)
library(gplots)
library(viridis)
library(pbapply)
library(pbmcapply)
library(ggrepel)
library(qvalue)
library(rtracklayer)
library(Gviz)
library(biomaRt)
library(cowplot)
library(stringr)

setwd("/scratch/groups/rmccoy22/syan11/sv_selection/ohana")
try(detach("package:plyr", unload = TRUE), silent = TRUE)

############################################################################

#system(cmd = "for i in {1..22}; do echo ${i}; zcat ~/work/syan11/sv_selection/ohana/k8/1KGP_SNP/vcfs/chr${i}.vcf.gz | grep -v ^# | cut -f8 | tr '=' '\t' | tr ';' '\t' | cut -f2,4 > ~/work/rmccoy22/sv_selection/pval/chr${i}_snp_af.txt ; done")

### Filter selscan results

pbs_results <- fread("pbs_results.txt")
gt_callrate <- fread("genotyping_callrate.txt")
hwe_pvals <- fread("HWE_pvals.txt")
all_pops_ld <- fread("/work-zfs/rmccoy22/syan11/sv_selection/snp_ld/all_pops_ld_20210114.txt") %>%
  .[, c("SV", "R2", "pop")]
svs <- fread("paragraph_sv_lengths.vcf")
afs <- fread(cmd = "zcat /scratch/groups/rmccoy22/rmccoy22/sv_selection/eichlerSVs_af_allpops_unfolded.frq.strat.gz")
common_svs <- unique(afs[MAF > 0.05 & MAF < 0.95]$SNP) # common in any 1KGP population

sv_freq <- afs %>%
  group_by(., SNP) %>%
  summarize(., af = sum(MAC) / sum(NCHROBS)) %>%
  as.data.table()
sv_freq[, maf := as.numeric(NA)]
sv_freq[af < 0.5, maf := af]
sv_freq[af >= 0.5, maf := 1 - af]

read_snp <- function(pop_number, chrom) {
  
  snp_freq <- fread(paste0("/scratch/groups/rmccoy22/rmccoy22/sv_selection/pval/chr", chrom, "_snp_af.txt"), header = FALSE) %>%
    setnames(., c("ac", "an")) %>%
    .[, af := ac / an]
  snp_freq[, maf := as.numeric(NA)]
  snp_freq[af < 0.5, maf := af]
  snp_freq[af >= 0.5, maf := 1 - af]
  
  selscan_snp <- fread(paste0("/scratch/groups/rmccoy22/syan11/sv_selection/ohana/k8/1KGP_SNP/selscan_rerun/chr",
                              chrom, "_selscan_k_p", pop_number, ".out"))
  snp_map <- fread(paste0("/scratch/groups/rmccoy22/syan11/sv_selection/ohana/k8/1KGP_SNP/vcfs/chr",
                          chrom, ".map"), header = FALSE) %>%
    setnames(., c("chr", "ID", "drop", "pos")) %>%
    .[, -3, with = FALSE]
  selscan_snp <- cbind(snp_map, snp_freq, selscan_snp)
  selscan_snp <- selscan_snp[`global-lle` > -1000] %>%
    setorder(., -`lle-ratio`)
  
  return(selscan_snp)
}

match_sv_snp_af <- function(selscan_snp, selscan_filt, multiple, maf_select) {
  n_selection <- nrow(selscan_filt[maf_bin == maf_select]) * multiple
  
  snp_selection <- selscan_snp[maf_bin == maf_select] %>%
    .[sample(1:.N, n_selection)]
  
  return(snp_selection)  
}

get_lle_perc <- function(freq_matched_snps, selscan_filt, row_index) {
  row_subset <- selscan_filt[row_index,]
  perc <- mean(freq_matched_snps$`lle-ratio` < row_subset$`lle_ratio`)
  return(perc)
}

filter_selscan <- function(p, common_svs_in, sv_freq_in, match_snp_p = FALSE) {
  selscan <- fread(paste0("k8/selscan/selscan_50_k8_p", p, ".out"))
  
  # restrict to autosomal SVs
  selscan_svs <- cbind(svs, selscan) %>%
    .[`#CHROM` %in% paste0("chr", 1:22)]
  
  # filter out SVs that violate HWE or have low callrate
  selscan_filt <- semi_join(selscan_svs,
                            gt_callrate[call_rate >= 0.5],
                            by = c("ID" = "SV")) %>%
    semi_join(., hwe_pvals[non_hwe_counts <= 13],
              by = c("ID" = "SV")) %>%
    as.data.table() %>%
    .[, target_pop := p]
  
  # filter out SVs that are not locally common
  selscan_filt <- selscan_filt[ID %in% common_svs_in]
  
  # filter out SVs with extreme global LLE
  selscan_filt <- selscan_filt[`global-lle` > -1000]
  
  setorder(selscan_filt, -`lle-ratio`)
  setnames(selscan_filt, "lle-ratio", "lle_ratio")
  selscan_filt[, rank := rank(-lle_ratio)]
  
  setnames(sv_freq_in, "SNP", "ID", skip_absent = TRUE)
  selscan_filt <- merge(selscan_filt, sv_freq_in, by = "ID")
  
  if (match_snp_p == TRUE) {
    ### compute percentile by comparing to SNPs and small indels
    
    #selscan_snp <- rbindlist(pblapply(as.character(c(1:22)),
    #                                  function(x) read_snp(p + 1, x)))
    selscan_snp <- read_snp(p + 1, "1")
    
    # match frequency distributions
    selscan_snp[, maf_bin := round(maf, 2)] %>%
      setorder(., maf_bin)
    selscan_filt[, maf_bin := round(maf, 2)] %>%
      setorder(., maf_bin)
    multiple <- floor(min(table(selscan_snp$maf_bin) / table(selscan_filt$maf_bin)))
    print(multiple)
    freq_matched_snps <- rbindlist(lapply(unique(selscan_filt$maf_bin),
                                          function(x) match_sv_snp_af(selscan_snp, selscan_filt, multiple, x)))
    
    snp_percentile <- pbmclapply(1:nrow(selscan_filt),
                       function(x) get_lle_perc(freq_matched_snps, selscan_filt, x),
                       mc.cores = getOption("mc.cores", 48)) %>%
      unlist()
    
    selscan_filt[, snp_perc := snp_percentile]
  }
  
  return(selscan_filt)
}

# snp-SV matching requires random sampling of SNPs, so set random seed
set.seed(1)
selscan_res <- rbindlist(pblapply(0:7, function(x) filter_selscan(x, common_svs, sv_freq, match_snp_p = TRUE))) %>%
  setorder(., -lle_ratio)

selscan_res[, p_nominal := pchisq(lle_ratio, df = 1, lower.tail = FALSE)]
selscan_res[, p_adj := p.adjust(p_nominal, method = "bonferroni")]

selscan_res[, snp_perc_99.9 := snp_perc > 0.999]
table(selscan_res$snp_perc_99.9)

############################################################################

setwd("/scratch/groups/rmccoy22/syan11/sv_selection/ohana")

### annotations: LD with 1KGP SNPs, genes from gencode

# add superpop ID columns to LD table
afr <- c("YRI", "LWK", "GWD", "MSL", "ESN", "ASW", "ACB")
amr <- c("MXL", "PUR", "CLM", "PEL")
eas <- c("CHB", "JPT", "CHS", "CDX", "KHV")
eur <- c("CEU", "TSI", "FIN", "GBR", "IBS")
sas <- c("GIH", "PJL", "BEB", "STU", "ITU")
all_pops_list <- list(afr, amr, eas, eur, sas) # all populations, grouped by superpop
all_superpops <- c("AFR", "AMR", "EAS", "EUR", "SAS")

get_superpops <- rbind(
  data.table(pop = afr, superpop = "AFR"),
  data.table(pop = amr, superpop = "AMR"),
  data.table(pop = eas, superpop = "EAS"),
  data.table(pop = eur, superpop = "EUR"),
  data.table(pop = sas, superpop = "SAS"))

### LD with SNPs
all_pops_ld <- merge(all_pops_ld, get_superpops, by = "pop")

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

exons <- gencode[type == "exon"]
exons_gr <- makeGRangesFromDataFrame(exons, 
                                     seqnames.field = "seqid", 
                                     start.field = "start", 
                                     end.field = "end",
                                     keep.extra.columns = TRUE)

# find overlaps between SVs and gene annotations
olaps <- findOverlaps(vcf_gr, genes_gr)
olaps_exons <- findOverlaps(vcf_gr, exons_gr)

genic_svs <- unique(vcf_gr[queryHits(olaps)]$ID)
exonic_svs <- unique(vcf_gr[queryHits(olaps_exons)]$ID)

length(unique(selscan_res[snp_perc_99.9 == TRUE & ID %in% genic_svs]$ID))
length(unique(selscan_res[snp_perc_99.9 == TRUE & ID %in% exonic_svs]$ID))

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
  dt_subset <- dt[target_pop == i]
  # get max LD of SV with a SNP in superpop of interest
  ld_subset <- all_pops_ld[superpop == superpop_group,] %>%
    setnames(., "R2", "max_R2_in_superpop") %>%
    # remove rows with duplicate SV names, keeping only the first copy that appears
    # (this is the copy with the highest r^2 value in that population)
    .[!duplicated(SV), ]
  setnames(ld_subset, "SV", "ID")
  
  # merge selscan results and LD data
  dt_subset <- merge(dt_subset, ld_subset,
                     by = "ID", all.x = TRUE)
  # merge selscan results and gene annotations
  dt_subset <- merge(dt_subset, sv_gene_overlaps,
                     by = "ID", all.x = TRUE)
  # add column specifying ancestry component number
  dt_subset$ancestry_component <- paste0("Ancestry component ", i + 1)
  # reorder by highest LLR
  setorder(dt_subset, -lle_ratio)
  
  return(dt_subset)
}

ancestry_to_superpop <- c("AMR","EAS","EUR","AMR","SAS","EAS","AFR","EUR")
selscan_res2 <- rbindlist(lapply(0:7, function(i) ld_genes_annot(selscan_res,
                                                                 ancestry_to_superpop[i + 1],
                                                                 i)))

selscan_res2[rank > 15, genes_string := NA]
length(unique(selscan_res2[snp_perc_99.9 == TRUE & max_R2_in_superpop > 0.8]$ID))

####

# overlap with eQTLs

read_qtl <- function(basedir, chunk) {
  message(chunk)
  dt <- fread(paste0(basedir, "qtls", chunk, ".out.permcovnew"), header = FALSE) %>%
    setnames(., c("gene_id", "nvars", "shape1", "shape2", "dummy", "sv_id", "dist", "nom_pval", "slope", "direct_perm_p", "beta_perm_p")) %>%
    .[, chunk := chunk]
  return(dt)
}

qtl <- rbindlist(lapply(1:22, function(x) read_qtl("/scratch/groups/rmccoy22/dnair4/QTLS_permutation_newcov/", x)))
qtl[, qval := qvalue(beta_perm_p)$qvalues]

selscan_res2[, is_eqtl := ID %in% unique(qtl[qval < 0.1]$sv_id)]

length(unique(selscan_res2[snp_perc_99.9 == TRUE & is_eqtl == TRUE]$ID))
qtl[sv_id %in% unique(selscan_res2[snp_perc_99.9 == TRUE & is_eqtl == TRUE]$ID) & qval < 0.1]

############################################################################

### plot of SV length vs. LLR

selscan_all <- selscan_res2

# no R^2 coloring for SVs that are not significant
selscan_all[p_adj > 0.01, max_R2_in_superpop := NA]

# no gene annotation for SVs that are not significant
selscan_all[is.na(genes_string), genes_string := ""]

save(selscan_all, file = "~/selscan_all.RData")
load("~/selscan_all.RData")

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

pdf(file = "~/Dropbox/papers/2020_sv/20201229_p2_R2.pdf", width = 20, height = 3)
p <- ggplot() +
  facet_wrap(~ ancestry_component, scales = "free", nrow = 1) +
  scale_x_log10() +
  theme_bw() +
  xlab("\nSV length (bp)") + ylab("Log-likelihood ratio statistic\n") +
  geom_point_rast(data = selscan_all, aes(x = abs(INFO), y = lle_ratio, color = max_R2_in_superpop), size = 0.5) +
  geom_text_repel(data = selscan_all[rank <= 15], 
                  aes(x = abs(INFO), y = lle_ratio, label = genes_string), 
                  size = 3, fontface = "italic", color = "black", force = 2) +
  scale_color_viridis_c(option = "plasma", name = "Maximum LD (r2)") +
  theme(legend.position = "none", panel.grid = element_blank()) +
  theme(panel.spacing = unit(1, "lines"))

g <- ggplot_gtable(ggplot_build(p))

strips <- which(grepl('strip-', g$layout$name))

pal <- gg_color_hue(8)

for (i in seq_along(strips)) {
  k <- which(grepl('rect', g$grobs[[strips[i]]]$grobs[[1]]$childrenOrder))
  g$grobs[[strips[i]]]$grobs[[1]]$children[[k]]$gp$fill <- pal[i]
}

plot(g)
dev.off()

###

table1 <- selscan_all[, c("ID", "INFO", "lle_ratio", "p_adj", "target_pop", "rank", "genes_string")] %>%
  setorder(., target_pop, rank) %>%
  .[, target_pop := target_pop + 1] %>%
  .[rank < 4] %>%
  .[, rank := NULL] %>%
  setnames(., c("SV ID", "SV length", "LRS", "Adj. p-value", "Ancestry component", "Affected gene(s)"))

formattable(table1, digits = 3, format = "markdown")
