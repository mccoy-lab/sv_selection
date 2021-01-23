library(data.table)
library(tidyverse)
library(GenomicRanges)
library(rtracklayer)
library(qvalue)
library(qqman)
library(ggrepel)
library(cowplot)

read_qtl <- function(basedir, chunk) {
  message(chunk)
  dt <- fread(paste0(basedir, "qtls", chunk, ".out.permcovnew"), header = FALSE) %>%
    setnames(., c("gene_id", "nvars", "shape1", "shape2", "dummy", "sv_id", "dist", "nom_pval", "slope", "direct_perm_p", "beta_perm_p")) %>%
    .[, chunk := chunk]
  return(dt)
}

qtl <- rbindlist(lapply(1:22, function(x) read_qtl("~/Downloads/", x)))


###

gt_callrate <- fread("~/Downloads/genotyping_callrate.txt")
hwe_pvals <- fread("~/Downloads/HWE_pvals.txt")

exclude <- unique(c(gt_callrate[call_rate < 0.5]$SV,
                    hwe_pvals[non_hwe_counts > 13]$SV))

qtl <- qtl[!(sv_id %in% exclude) & !is.na(slope)]
qtl[, sv_type := toupper(sapply(strsplit(sv_id, "_"), "[[", 3))]
qtl[, qval := qvalue(beta_perm_p)$qvalues]

###

gencode <- data.table(readGFF("ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_26/gencode.v26.annotation.gtf.gz"))
genes <- gencode[type == "gene"][, c("seqid", "start", "end", "gene_id", "gene_name", "gene_type")]
qtl <- merge(qtl, genes, "gene_id")

setorder(qtl, beta_perm_p)
qtl[, rank := .I]
qtl[, expected_p := rank / .N]

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
pal <- gg_color_hue(1)

p1 <- ggplot(data = qtl, aes(x = -log10(expected_p), y = -log10(beta_perm_p), color = qval < 0.1)) +
  geom_point() +
  theme_bw() +
  geom_abline(slope = 1) +
  theme(panel.grid = element_blank(),
        legend.position = "none") +
  scale_color_manual(values = c("grey35", "#9E129D")) +
  labs(y = expression(observed-log[10](p)),
       x = expression(expected-log[10](p)))

p2 <- ggplot() +
  geom_point(data = qtl, aes(x = slope, y = -log10(beta_perm_p), color = qval < 0.1)) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = "none") +
  scale_color_manual(values = c("grey35", "#9E129D")) +
  labs(y = expression(-log[10](p)),
       x = expression(beta)) +
  geom_text_repel(data = qtl[beta_perm_p < 1e-30], 
                  aes(x = slope, y = -log10(beta_perm_p), label = gene_name),
                  fontface = "italic", size = 3)

p3 <- ggplot(data = qtl[qval < 0.1], aes(x = dist)) +
  geom_histogram() +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  xlab("SV distance from TSS (bp)") +
  ylab("Number of eQTLs")

plot_grid(p1, p3, p2, nrow = 1, 
          labels = c("A.", "B.", "C."),
          rel_widths = c(1, 1.05, 1.1))
