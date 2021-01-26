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

read_snp <- function(pop_number, chrom, minimal = FALSE) {
  
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
  
  if (minimal == TRUE) {
    selscan_snp <- selscan_snp[, c("ID", "af", "maf", "lle-ratio")]
  }
  
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
    
    # get chromosome 1 SNP data from same ancestry component
    selscan_snp <- read_snp(p + 1, "1", minimal = FALSE)
    
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

length(unique(selscan_res[snp_perc_99.9 == TRUE]$ID))

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


############################################################################

# plot locus-zoom style local manhattan plot

get_locuszoom <- function(sv_id, pop_number, selscan_results_in, window_size) {
  
  selscan_subset <- selscan_results_in[ID == sv_id & target_pop == pop_number] %>%
    .[, 1:16] %>%
    .[, is_sv := TRUE] %>%
    setnames(., "#CHROM", "chr") %>%
    .[, chr := gsub("chr", "", chr)] %>%
    setnames(., "POS", "pos")
  chrom <- selscan_subset$chr
  position <- selscan_subset$pos
  selscan_snp <- read_snp(pop_number + 1, chrom) %>%
    .[pos > (position - window_size / 2) & pos < (position + window_size / 2)] %>%
    .[, is_sv := FALSE] %>%
    .[, -c("ac", "an", "af", "maf")] %>%
    .[, INFO := 1] %>%
    setnames(., "lle-ratio", "lle_ratio")
  
  keep_col <- colnames(selscan_snp)
  
  selscan_local <- rbind(selscan_subset[, ..keep_col], selscan_snp) %>%
    .[, locus_title := paste0("Ancestry component ", pop_number + 1)]
  
  return(selscan_local)
}


locus_0 <- get_locuszoom("27407_HG02106_ins", 0, selscan_res, 4e7)
locus_1 <- get_locuszoom("22237_HG02059_ins", 1, selscan_res, 1e6)
locus_2 <- get_locuszoom("32021_HG00268_ins", 2, selscan_res, 2e6)
locus_3 <- get_locuszoom("22065_HG02106_del", 3, selscan_res, 3e6)
locus_4 <- get_locuszoom("21859_NA19240_ins", 4, selscan_res, 3e6)
locus_5 <- get_locuszoom("22237_HG02059_ins", 5, selscan_res, 1e6)
locus_6 <- get_locuszoom("10085_HG00268_del", 6, selscan_res, 1e6)
locus_7 <- get_locuszoom("18075_HG00268_del", 7, selscan_res, 1e6)

locuszoom <- rbind(locus_0, locus_1, locus_2, locus_3, locus_4, locus_5, locus_6, locus_7)
locuszoom <- rbind(locuszoom[lle_ratio == 0][sample(1:nrow(locuszoom[lle_ratio == 0]), 1e5)],
                   locuszoom[lle_ratio != 0]) # sparse the points at 0 for plotting efficiency

ggplot(data = locuszoom[is_sv == FALSE], aes(x = pos / 1000000, y = lle_ratio)) + 
  scale_color_manual(values = "black") +
  geom_point(size = 0.5) +
  geom_point(data = locuszoom[is_sv == TRUE], color = "red", size = 2) +
  geom_text_repel(data = locuszoom[is_sv == TRUE], aes(label = ID), 
                  color = "black", size = 2) +
  ylab("Likelihood ratio statistic") +
  xlab("Position (Mbp)") +
  theme_bw() +
  theme(legend.position = "none", 
        panel.grid = element_blank()) +
  facet_wrap(~ locus_title, scales = "free")

############################################################################

# local investigation of the IGHG4 region

ighg4 <- get_locuszoom("22237_HG02059_ins", 1, selscan_res, 1e6)

check_neand_vcf <- function(locus, row_index) {
  tryCatch({
    snp <- locus[row_index]
    ref <- strsplit(snp$ID, "_")[[1]][3]
    alt <- strsplit(snp$ID, "_")[[1]][4]
    if (nchar(ref) != 1 || nchar(alt) != 1) {
      return(NA)
    }
    cmd <- paste("/home-net/home-4/rmccoy22@jhu.edu/code/htslib-1.11/tabix",
                 paste0("/scratch/groups/rmccoy22/rmccoy22/sv_selection/neand/AltaiNea.hg19_1000g.", snp$chr, ".mod.vcf.gz"),
                 paste0(snp$chr, ":", snp$hg19_pos, "-", snp$hg19_pos))
    test <- strsplit(system(command = cmd, intern = TRUE), "\t")
    if (test[[1]][5] == alt & (grepl("^1/1", test[[1]][10]) | (grepl("^0/1", test[[1]][10])))) {
      return(TRUE)
    } else {
      return(FALSE)
    }
  }, error = function(error_condition) {
    return(FALSE)
  })
}

check_archaic_bam <- function(locus, row_index, archaic) {
  tryCatch({
    snp <- locus[row_index]
    ref <- strsplit(snp$ID, "_")[[1]][3]
    alt <- strsplit(snp$ID, "_")[[1]][4]
    chrom <- unique(locus$chr)
    if (nchar(ref) != 1 || nchar(alt) != 1) {
      return(NA)
    }
    if (archaic == "altai") {
      cmd <- paste(paste0("/home-net/home-4/rmccoy22@jhu.edu/code/samtools-1.11/samtools mpileup AltaiNea.hg19_1000g.", chrom, ".dq.bam -r"),
                   paste0(snp$chr, ":", snp$hg19_pos, "-", snp$hg19_pos))
    } else if (archaic == "vindija") {
      cmd <- paste(paste0("/home-net/home-4/rmccoy22@jhu.edu/code/samtools-1.11/samtools mpileup Vi33.19.chr", chrom, ".indel_realn.bam -r"),
                   paste0(snp$chr, ":", snp$hg19_pos, "-", snp$hg19_pos))
    } else if (archaic == "chagyrskaya") {
      cmd <- paste(paste0("/home-net/home-4/rmccoy22@jhu.edu/code/samtools-1.11/samtools mpileup chr", chrom, ".rh.bam -r"),
                   paste0(snp$chr, ":", snp$hg19_pos, "-", snp$hg19_pos))
    } else if (archaic == "denisova") {
      cmd <- paste("/home-net/home-4/rmccoy22@jhu.edu/code/samtools-1.11/samtools mpileup T_hg19_1000g.bam -r",
                   paste0(snp$chr, ":", snp$hg19_pos, "-", snp$hg19_pos))
    } else {
      return(NA)
    }
    query_results <- strsplit(system(command = cmd, intern = TRUE), "\t")
    if (str_count(toupper(query_results[[1]][5]), alt) > 1) {
      return(TRUE)
    } else {
      return(FALSE)
    }
  }, error = function(error_condition) {
    return(NA)
  })
}

ighg4[, chrom := paste0("chr", chr)]
ighg4[, index := .I]
hg38_coords <- makeGRangesFromDataFrame(ighg4, 
                                        seqnames.field = "chrom", 
                                        start.field = "pos", 
                                        end.field = "pos", 
                                        ignore.strand = TRUE)

chain <- import.chain("/scratch/groups/rmccoy22/rmccoy22/sv_selection/neand/hg38ToHg19.over.chain")
hg19_coords <- liftOver(hg38_coords, chain) %>%
  as.data.table() %>%
  setnames(., "group", "index") %>%
  setnames(., "start", "hg19_pos")

ighg4 <- merge(ighg4, hg19_coords[, c("index", "hg19_pos")], by = "index")

# download bam indices
setwd("/scratch/groups/rmccoy22/rmccoy22/sv_selection/neand/")

# get Altai BAM
url <- "http://cdna.eva.mpg.de/neandertal/altai/AltaiNeandertal/bam/AltaiNea.hg19_1000g.14.dq.bam"
download.file(url, destfile = basename(url))
url <- "http://cdna.eva.mpg.de/neandertal/altai/AltaiNeandertal/bam/AltaiNea.hg19_1000g.14.dq.bam.bai"
download.file(url, destfile = basename(url))

# get Vindija BAM
url <- "http://cdna.eva.mpg.de/neandertal/Vindija/bam/Vi33.19.chr14.indel_realn.bam"
download.file(url, destfile = basename(url))
url <- "http://cdna.eva.mpg.de/neandertal/Vindija/bam/Vi33.19.chr14.indel_realn.bam.bai"
download.file(url, destfile = basename(url))

# get Chagyrskaya BAM
url <- "http://cdna.eva.mpg.de/neandertal/Chagyrskaya/bam/chr14.rh.bam"
download.file(url, destfile = basename(url))
url <- "http://cdna.eva.mpg.de/neandertal/Chagyrskaya/bam/chr14.rh.bam.bai"
download.file(url, destfile = basename(url))

# get Denisova BAM
url <- "http://cdna.eva.mpg.de/denisova/alignments/T_hg19_1000g.bam"
download.file(url, destfile = basename(url))
url <- "http://cdna.eva.mpg.de/denisova/alignments/T_hg19_1000g.bam.bai"
download.file(url, destfile = basename(url))

# compare to each of the archaic genome alignments
altai_match <- unlist(pbmclapply(1:nrow(ighg4), function(x) check_archaic_bam(ighg4, x, archaic = "altai"), 
                                 mc.cores = getOption("mc.cores", 48)))

vndja_match <- unlist(pbmclapply(1:nrow(ighg4), function(x) check_archaic_bam(ighg4, x, archaic = "vindija"), 
                                 mc.cores = getOption("mc.cores", 48)))

chgyr_match <- unlist(pbmclapply(1:nrow(ighg4), function(x) check_archaic_bam(ighg4, x, archaic = "chagyrskaya"), 
                                 mc.cores = getOption("mc.cores", 48)))

denis_match <- unlist(pbmclapply(1:nrow(ighg4), function(x) check_archaic_bam(ighg4, x, archaic = "denisova"), 
                                 mc.cores = getOption("mc.cores", 48)))

ighg4[, match_altai := altai_match]
ighg4[, match_vndja := vndja_match]
ighg4[, match_denis := denis_match]
ighg4[, match_chgyr := chgyr_match]

# get SPrime calls
url <- "https://data.mendeley.com/public-files/datasets/y7hyt83vxr/files/f01566ae-0a9b-4847-b312-1af059fae3d1/file_downloaded"
download.file(url, destfile = basename(url))
system(paste("tar -zxvf", basename(url)))
sprime <- fread("mendeley_data/CDX.chr14.ND_match")
ighg4[, is_sprime_sig := hg19_pos %in% sprime$POS]

min(ighg4[lle_ratio > 450]$pos)
max(ighg4[lle_ratio > 450]$pos)

############################################################################

# generate haplostrips

# /home-net/home-4/rmccoy22@jhu.edu/code/htslib-1.11/tabix -h http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20181203_biallelic_SNV/ALL.chr14.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.gz chr14:105498101-105822843 > ighg4_1kgp.vcf

haplostrips <- ighg4[lle_ratio > 450]

kgp <- fread("ighg4_1kgp.vcf") %>%
  .[POS %in% haplostrips$pos] %>%
  .[, ID := paste0(`#CHROM`, "_", POS, "_", REF, "_", ALT)] %>%
  .[, c(3, 10:(ncol(.) - 1)), with = FALSE] %>%
  pivot_longer(!ID, names_to = "sample_id", values_to = "gt") %>%
  as.data.table() %>%
  .[, c("h1", "h2") := tstrsplit(gt, "|", fixed=TRUE)] %>%
  .[, gt := NULL] %>%
  pivot_longer(!c(ID, sample_id), names_to = "haplotype", values_to = "gt") %>%
  as.data.table() %>%
  rbind(., data.table(ID = haplostrips$ID, sample_id = "altai", haplotype = NA, gt = as.numeric(haplostrips$match_altai))) %>%
  rbind(., data.table(ID = haplostrips$ID, sample_id = "vindija", haplotype = NA, gt = as.numeric(haplostrips$match_vndja))) %>%
  rbind(., data.table(ID = haplostrips$ID, sample_id = "chagyrskaya", haplotype = NA, gt = as.numeric(haplostrips$match_chgyr))) %>%
  rbind(., data.table(ID = haplostrips$ID, sample_id = "denisova", haplotype = NA, gt = as.numeric(haplostrips$match_denis)))

sample_manifest <- fread("igsr_samples.tsv") %>%
  .[, c(1, 4, 6)] %>%
  setnames(., c("sample_id", "pop", "superpop"))

kgp <- merge(kgp, sample_manifest, by = "sample_id", all.x = TRUE, allow.cartesian = TRUE)
kgp[sample_id %in% c("altai", "vindija", "chagyrskaya", "denisova"), pop := "archaic"]
kgp[sample_id %in% c("altai", "vindija", "chagyrskaya", "denisova"), superpop := "archaic"]

kgp[, pos := sapply(strsplit(kgp$ID, "_"), "[[", 2)]
setorder(kgp, pos)

kgp$ID <- factor(kgp$ID, levels = kgp[!duplicated(pos)]$ID)

kgp_keep <- kgp[pop %in% c("archaic", "CDX", "KHV", "CHB", "JPT", "CEU", "YRI") & !is.na(gt)] %>%
  .[, hap_id := paste(sample_id, haplotype, sep = "_")]
kgp_keep$pop <- factor(kgp_keep$pop, levels =  c("archaic", "CDX", "KHV", "CHB", "JPT", "CEU", "YRI"))

sample_order_dt <- kgp_keep[, -"pos"] %>%
  .[, -c("sample_id", "haplotype", "pop", "superpop")] %>%
  pivot_wider(., names_from = hap_id, values_from = gt, values_fn = list(gt = unique)) %>%
  as.data.table() %>%
  .[, -"ID"] %>%
  t(.)

orig_order <- rownames(sample_order_dt)

sample_order <- sample_order_dt %>%
  dist(., method = "binary") %>%
  hclust(.)

new_order <- orig_order[sample_order$order]

kgp_keep$hap_id <- factor(kgp_keep$hap_id, levels = new_order)

keep_samples <- group_by(kgp_keep[pop != "archaic"], hap_id) %>%
  summarize(., pop = unique(pop)) %>%
  group_by(., pop) %>%
  summarize(., keep_list = list(sample(hap_id, 30))) %>%
  as.data.table() %>%
  .$keep_list %>%
  unlist()

keep_samples <- c(as.character(keep_samples), "altai_NA", "chagyrskaya_NA", "denisova_NA", "vindija_NA")

# haplostrips
ggplot(data = kgp_keep[hap_id %in% keep_samples], aes(x = ID, y = sample_id, fill = gt)) + 
  geom_tile() +
  scale_fill_manual(values = c("white", "black")) +
  theme(axis.text = element_blank(), 
        legend.position = "none", 
        panel.spacing = unit(0, "lines"),
        panel.grid = element_blank(),
        panel.background = element_blank()) +
  facet_grid(pop ~ ., scales = "free") +
  xlab("SNP") +
  ylab("Sample") +
  NULL

# local LRS plot
ggplot(data = ighg4[is_sv == FALSE & !is.na(match_chgyr)], aes(x = pos / 1000, y = lle_ratio, color = match_chgyr)) + 
  scale_color_manual(values = c("grey55", "purple")) +
  geom_point(size = 0.5) +
  geom_point(data = ighg4[is_sv == TRUE & is.na(match_chgyr)], color = "red", size = 2) +
  ylab("Likelihood ratio statistic") +
  xlab("Position (Kbp)") +
  theme_bw() +
  theme(legend.position = "none", 
        panel.grid = element_blank()) +
  geom_hline(yintercept = 450, lty = "dashed", alpha = 0.2)                        
