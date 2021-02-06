library(data.table)
library(dplyr)

### USAGE: Rscript add_NA18498.R <VCF>
### Adds genotypes for the NA18498 sample, which is missing from the 1KGP
### GRCh38 variant calls. These genotypes are copied from the NA18499 sample.

############################################################################


# accept VCF filepath as arguments
args <- commandArgs(trailingOnly = TRUE)
# check if VCF was provided
if ( length(args) < 1 ) {
  stop("Input VCF is required",
       call. = FALSE)
}

vcfpath <- args[1]

# read in VCF
vcf <- fread(cmd = paste0("grep -v '##' ", vcfpath)) %>%
  # create "NA18498" column and assign it genotypes from NA18499
  .[, NA18498 := NA18499] %>%

# output modified VCF (without a header) to file called `NA18498.vcf`
fwrite(vcf, "NA18498.vcf",
       sep = "\t", row.names = FALSE, col.names = TRUE)