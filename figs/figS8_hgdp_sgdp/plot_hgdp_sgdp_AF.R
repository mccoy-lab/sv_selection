library(data.table)
library(dplyr)
library(ggplot2)
library(patchwork)

### Plot global AFs of a SNP in HGDP and SGDP.

#########################################################################

### Convert HGDP metadata into a `within` file for calculating AF with plink.
# metadata downloaded from: https://www.internationalgenome.org/data-portal/data-collection/hgdp

# HGDP
hgdp_meta <- fread("igsr_hgdp.tsv")
within_file <- hgdp_meta[, c(1, 1, 8)] %>%
  .[`Population elastic ID` == "PapuanSGDP,PapuanHighlandsHGDP", `Population elastic ID` := "PapuanHGDP"] %>%
  .[`Population elastic ID` == "PapuanSGDP,PapuanSepikHGDP", `Population elastic ID` := "PapuanHGDP"] %>%
  .[`Population elastic ID` == "PapuanSepikHGDP", `Population elastic ID` := "PapuanHGDP"]
fwrite(within_file, "hgdp_within.txt",
       sep = "\t", col.names = FALSE)
# calculate AF with ./calculage_hdgp.af.sh


#########################################################################

### Plot AFs by population.

### HGDP
# use metadata to determine how to group populations
# also clean up population names
pop_meta <- unique(hgdp_meta[,c(7,8)]) %>%
  .[`Population elastic ID` == "PapuanSGDP,PapuanHighlandsHGDP", `Population elastic ID` := "PapuanHGDP"] %>%
  .[`Population elastic ID` == "PapuanSGDP,PapuanSepikHGDP", `Population elastic ID` := "PapuanHGDP"] %>%
  .[`Superpopulation name` == "Oceania (SGDP),Oceania (HGDP)", `Superpopulation name` := "Oceania"] %>%
  .[`Population elastic ID` == "PapuanSepikHGDP", `Population elastic ID` := "PapuanHGDP"]
setnames(pop_meta, c("Superpopulation name", "Population elastic ID"), c("superpop", "pop"))
pop_meta$superpop <- gsub(" \\(HGDP\\)", "", pop_meta$superpop)
pop_meta <- unique(pop_meta) # remove duplicate populations

# read in AFs, select only SNP of interest
hgdp <- fread("hgdp_chr14.frq.strat") %>%
  .[SNP == "rs150526114"] %>%
  # calculate alt AF
  .[, alt_AF := 1 - MAF]
# merge with metadata
hgdp <- merge(hgdp, pop_meta,
              by.x = "CLST", by.y = "pop")
hgdp$CLST <- gsub("HGDP", "", hgdp$CLST) # clean up names
# make separate alt AF and ref AF rows for barplotting
hgdp <- gather(hgdp, ident, AF, MAF, alt_AF)

# plot
plot_hgdp <- function(pop_vec) {
  subset <- hgdp %>%
    setDT() %>%
    .[superpop %in% pop_vec]
  subset$superpop <- factor(subset$superpop,
                            levels = pop_vec)
  
  p <- ggplot(subset,
              aes(fill = ident, y = AF, x = CLST)) + 
    geom_bar(position = "stack", stat = "identity") +
    scale_fill_manual("Allele:", values = c("MAF" = "#9E139D", "alt_AF" = "#FFC01D"), labels = c("REF", "ALT")) +
    # xlab("Population") + ylab("Allele frequency") +
    theme_bw() +
    theme(panel.border = element_blank(), panel.grid = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1),
          # legend.position = "none",
          plot.title = element_text(hjust = 0.5),
          axis.title.x = element_blank(), axis.title.y = element_blank()) +
    facet_grid(~ superpop, space = "free", scales = "free_x")
  
  return(p)
}
# generate plots
hgdp1 <- plot_hgdp(c("East Asia",  "Oceania", "Central South Asia")) +
  ggtitle("HGDP")
hgdp2 <- plot_hgdp(c("Africa", "Middle East", "Europe", "America"))

### SGDP
# read in AFs, select only SNP of interest
sgdp <- fread("sgdp_afs_population.txt") %>%
  .[Variant == "chr14_105613313_G_A"] %>%
  .[Region == "EastAsia", Region := "East Asia"] %>%
  .[Region == "SouthAsia", Region := "South Asia"] %>%
  .[Region == "CentralAsiaSiberia", Region := "Central Asia Siberia"] %>%
  .[Region == "WestEurasia", Region := "West Eurasia"]

# plot
plot_sgdp <- function(pop_vec) {
  subset <- sgdp[Region %in% pop_vec]
  subset$Region <- factor(subset$Region,
                          levels = pop_vec)
  
  p <- ggplot(subset,
              aes(fill = ident, y = AF, x = Population, order = ident)) + 
    geom_bar(position = position_fill(reverse = TRUE), stat = "identity") +
    scale_fill_manual("Allele:", values = c("REF" = "#FFC01D", "ALT" = "#9E139D")) +
    theme_bw() +
    theme(panel.border = element_blank(), panel.grid = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title = element_text(hjust = 0.5),
          axis.title.x = element_blank(), axis.title.y = element_blank()) +
    guides(fill = guide_legend(reverse = TRUE)) +
    facet_grid(~ Region, space = "free", scales = "free_x")
  
  return(p)
}
# generate plots
sgdp1 <- plot_sgdp(c("East Asia",  "Oceania", "Central Asia Siberia")) + 
  ggtitle("\nSGDP")
sgdp2 <- plot_sgdp(c("South Asia", "Africa"))
sgdp3 <- plot_sgdp(c("West Eurasia", "America"))

### combine plots using patchwork
p <- (hgdp1 / hgdp2 / sgdp1 / sgdp2 / sgdp3) +
  plot_layout(guides = "collect") & theme(legend.position = 'bottom')

ggsave("chr14_105613313_G_A.pdf", p, width = 7, height = 13)