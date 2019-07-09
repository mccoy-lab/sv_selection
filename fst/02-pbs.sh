library(data.table)
library(magrittr)
library(ggplot2)
library(dplyr)
library(pbapply)
library(qqman)

# script to read in pairwise FSTs and compute PBS
# written for previous project and needs to be generalized

pbs <- function(chrom, prefix) {
  rps_chb <- fread(paste(prefix, "chr", chrom, ".rps.chb.fst", sep = "")) %>%
    setnames(., c("chr", "pos", "rps.af", "chb.af", "fst.rps.chb"))
  rps_png <- fread(paste(prefix, "chr", chrom, ".rps.png.fst", sep = "")) %>%
    setnames(., c("chr", "pos", "rps.af", "png.af", "fst.rps.png"))
  png_chb <- fread(paste(prefix, "chr", chrom, ".png.chb.fst", sep = "")) %>%
    setnames(., c("chr", "pos", "png.af", "chb.af", "fst.png.chb"))
  fst <- merge(merge(rps_chb, rps_png, "pos"), png_chb, "pos")
  fst$pbs <- ((-log(1 - fst$fst.rps.png)) + (-log(1 - fst$fst.rps.chb)) - (-log(1 - fst$fst.png.chb))) / 2
  setorder(fst, -pbs)
  return(data.table(chr = fst$chr.x, pos = fst$pos, rps.af = fst$rps.af.x,
                    png.af = fst$png.af.x, chb.af = fst$chb.af.x, pbs = fst$pbs))
}

pbs_window <- function(chrom, prefix) {
  rps_chb <- fread(paste(prefix, "chr", chrom, ".rps.chb.windowed.weir.fst", sep = ""))
  setnames(rps_chb, "MEAN_FST", "fst.rps.chb")
  png_chb <- fread(paste(prefix, "chr", chrom, ".png.chb.windowed.weir.fst", sep = ""))
  setnames(png_chb, "MEAN_FST", "fst.png.chb")
  rps_png <- fread(paste(prefix, "chr", chrom, ".png.rps.windowed.weir.fst", sep = ""))
  setnames(rps_png, "MEAN_FST", "fst.rps.png")
  fst <- merge(merge(rps_chb, rps_png, "BIN_START"), png_chb, "BIN_START")
  fst[, pbs := (-log(1 - fst.rps.png)) + (-log(1 - fst.rps.chb)) - (-log(1 - fst.png.chb)) / 2]
}

pbs_results <- do.call(rbind, lapply(1:22, function(x) pbs(x, "/home/rajivm/home_vol2/flores.nobackup/FST_RPS_PNG_CHB/")))
pbs_results <- pbs_results[pbs > -Inf & pbs < Inf]
pbs_results <- pbs_results[rps.af != 0 & png.af != 0 & chb.af != 0] # require variant to be present in all populations
pbs_results[, mergeID := paste(chr, pos, sep = "_")]
setorder(pbs_results, chr, pos)

summarize_pbs <- function(pbs_window) {
  return(data.table(chrom = unique(pbs_window$chr), pos = median(pbs_window$pos), rps.af = mean(pbs_window$rps.af), 
                    png.af = mean(pbs_window$png.af), chb.af = mean(pbs_window$chb.af),
                    mean_pbs = mean(pbs_window$pbs), max_pbs = max(pbs_window$pbs)))
}

slide_window <- function(pbs_results, chrom, window_width, increment) {
  pbs_subset <- pbs_results[chr == chrom]
  n_snps <- nrow(pbs_subset)
  window_list <- seq(1, n_snps - window_width, increment)
  chrom_results <- do.call(rbind, pblapply(window_list, function(x) summarize_pbs(pbs_subset[x:(x + window_width),])))
}

pbs_window <- do.call(rbind, pblapply(1:22, function(x) slide_window(pbs_results, x, 20, 5)))

ihs <- fread("/net/akey/vol1/home/rcmccoy/vol2home/flores.nobackup/iHS/ihs.out") %>%
  setnames(., c("chr", "pos", "ihh0", "ihh1", "ihs", "ihs_std"))
	
ihs[, mergeID := paste(chr, pos, sep = "_")]
ihs[, chr := NULL]
ihs[, pos := NULL]

sel <- merge(pbs_results, ihs, "mergeID")

standardize_ihs <- function(ihs_table, af) {
	table_subset <- ihs_table[rps.af == af]
	table_subset$ihs_std <- scale(table_subset$ihs)
	return(table_subset)
}

sel_rescale <- do.call(rbind, lapply(unique(sel$rps.af), function(x) standardize_ihs(sel, x))) %>%
  setorder(., pbs)

sel_rescale[pbs < 0]$pbs <- 0
	
png("~/test.png", type = "cairo")
ggplot(data = sel_rescale, aes(x = pbs, y = abs(ihs_std), label = merge.id)) + 
  geom_text(size = 3) + 
	theme_bw() +
	xlab("PBS") +
	ylab("|iHS|")
dev.off()

write.table(sel_rescale, file = "~/vol2home/flores.nobackup/sel_rescale.txt", sep = "\t", row.names = F, col.names = T, quote = F)

sel_genes <- fread("/net/akey/vol1/home/rcmccoy/vol2home/flores.nobackup/sel_rescale.genes.bed") %>%
  setnames(., c("chr", "pos", "pos.x", "rps.af", "png.af", "chb.af", "pbs", "ihs_std", "gene.id", "gene.symbol"))
	
het_fp <- fread("~/vol2home/flores.nobackup/filters/png.HET_filter_AnyInd_sorted_22_coord.bed") %>%
	setnames(., c("chr", "pos", "pos.x", "gt"))
	
sel_genes[, merge.id := paste(chr, pos, sep = "_")]
het_fp[, merge.id := paste(chr, pos, sep = "_")]

sel_genes <- sel_genes[!(merge.id %in% het_fp$merge.id)]

png("~/test.png", type = "cairo")
ggplot(data = sel_genes, aes(x = pbs, y = abs(ihs_std), label = gene.symbol)) + 
	geom_text(size = 3) + 
  theme_bw() +
	xlab("PBS") +
	ylab("|iHS|")
dev.off()

pdf("~/test.pdf")
ggplot(data = sel_genes, aes(x = pbs, y = abs(ihs_std), label = gene.symbol)) + 
	geom_text(size = 3) + 
  theme_bw() +
	xlab("PBS") +
	ylab("|iHS|")
dev.off()

write.table(unique(sel_genes[pbs > 0.5170466 & abs(ihs_std) > 2.650268]$gene.id), file = "~/test.txt", quote = F, row.names = F, col.names = F)
write.table(unique(sel_genes$gene.id), file = "~/background.txt", quote = F, row.names = F, col.names = F)
