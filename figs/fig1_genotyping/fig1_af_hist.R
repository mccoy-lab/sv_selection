library(data.table)
library(tidyverse)
library(ggplot2)
library(sjPlot)
library(scales)

#########################################################################
### Read in AFs from PLINK and plot histograms by continental population.
#########################################################################

setwd("~/Downloads/")

gt_callrate_path <- "genotyping_callrate.txt"
hwe_path <- "HWE_pvals.txt"
af_path <- "eichlerSVs_af_allpops_unfolded.frq.strat.gz"

afs <- fread(af_path)
gt_callrate <- fread(gt_callrate_path)
hwe_pvals <- fread(hwe_path)

afs <- semi_join(afs,
                 gt_callrate[ call_rate >= 0.5 ],
                 by = c("SNP" = "SV")) %>%
  semi_join(., hwe_pvals[ non_hwe_counts <= 13 ],
            by = c("SNP" = "SV")) %>%
  as.data.table() %>%
  setnames("CLST", "pop")

# group SVs by continental population
afr <- c("YRI", "LWK", "GWD", "MSL", "ESN", "ASW", "ACB")
amr <- c("MXL", "PUR", "CLM", "PEL")
eas <- c("CHB", "JPT", "CHS", "CDX", "KHV")
eur <- c("CEU", "TSI", "FIN", "GBR", "IBS")
sas <- c("GIH", "PJL", "BEB", "STU", "ITU")

get_superpops <- rbind(
  data.table(pop = afr, superpop = "AFR"),
  data.table(pop = amr, superpop = "AMR"),
  data.table(pop = eas, superpop = "EAS"),
  data.table(pop = eur, superpop = "EUR"),
  data.table(pop = sas, superpop = "SAS"))

afs <- merge(afs, get_superpops, by = "pop")

superpop_afs <- group_by(afs, SNP, superpop) %>%
  summarize(sum_MAC = sum(MAC), sum_NCHROBS = sum(NCHROBS), maf = sum(MAC) / sum(NCHROBS)) %>%
  as.data.table()

global_afs <- group_by(afs, SNP) %>%
  summarize(sum_MAC = sum(MAC), sum_NCHROBS = sum(NCHROBS), maf = sum(MAC) / sum(NCHROBS)) %>%
  as.data.table()

ggplot(data = superpop_afs, 
       aes(x = maf, fill = superpop)) +
  geom_histogram() +
  facet_grid(. ~ superpop) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = "none") +
  xlab("Alternative Allele Frequency") +
  ylab("Number of SVs")

sv_metadata <- fread("paragraph_sv_lengths.vcf") %>%
  setnames(., "ID", "SNP")

global_afs <- merge(global_afs, sv_metadata, by = "SNP")
global_afs[, MAF := as.numeric(NA)]
global_afs[maf < 0.5, MAF := maf]
global_afs[maf >= 0.5, MAF := 1 - maf]

global_afs[, sv_type := sapply(strsplit(global_afs$SNP, "_"), "[", 3)]
global_afs[, sv_length := abs(INFO)]

ggplot(data = global_afs[MAF != 0],
       aes(y = abs(INFO), x = MAF, color = sv_type)) +
  geom_point(size = 0.1) +
  scale_y_log10() +
  scale_x_log10() +
  stat_smooth(method = "lm", formula = y ~ exp(x)) +
  facet_grid(. ~ sv_type)

global_afs[sv_type == "del", sv_type := "deletion"]
global_afs[sv_type == "dup", sv_type := "duplication"]
global_afs[sv_type == "ins", sv_type := "insertion"]

ggplot(data = global_afs[MAF != 0], 
       aes(x = MAF, fill = sv_type)) +
  geom_histogram() +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = "none") +
  xlab("Minor Allele Frequency") +
  ylab("Number of SVs")  +
  facet_grid(. ~ sv_type)

ggplot(data = global_afs[MAF != 0],
       aes(x = abs(INFO), y = MAF)) +
  geom_bin2d() +
  scale_y_log10() +
  scale_x_log10() +
  facet_grid(. ~ sv_type) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  scale_fill_viridis_c(name = "Number of SVs") +
  ylab("Minor allele frequency") +
  xlab("SV length (bp)")

cor.test(abs(global_afs[MAF != 0 & MAF != 1]$INFO),
         global_afs[MAF != 0 & MAF != 1]$MAF,
         method = "spearman")

lm(formula = log(MAF) ~ sv_type, data = global_afs[MAF != 0]) %>%
  summary()

vcf <- fread("cat ~/Downloads/origVariants_post50bpMerge_withSupp.vcf | grep -v '^##' | cut -f1-3,8")
vcf[, info_list := list(strsplit(INFO, ";"))]
sinfo <- lapply(vcf$info_list, tail, n = 2L)
supp <- as.numeric(gsub("SUPP=", "", unlist(lapply(sinfo, `[[`, 1))))
vcf[, supp := supp]

global_afs[, ID := gsub("_del", "", ID)]
global_afs[, ID := gsub("_ins", "", ID)]
global_afs[, ID := gsub("_dup", "", ID)]

setnames(global_afs, "SNP", "ID")
global_afs <- merge(global_afs, vcf[, c("ID", "supp")], by = "ID")

ggplot(data = global_afs, aes(x = factor(supp), y = maf)) +
  geom_boxplot(fill = "#85BAD0") +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  xlab("Number of LR-seq samples") +
  ylab("Paragraph 1KGP allele frequency")








