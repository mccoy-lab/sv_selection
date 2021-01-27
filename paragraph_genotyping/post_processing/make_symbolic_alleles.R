#!/usr/bin/env Rscript

library(dplyr)
library(data.table)
library(pbmcapply)

### USAGE: ./02_make_symbolic_alleles.R <input VCF> <manifest file>
###
### Replaces full SV sequence in ref/alt fields of VCF with symbolic representations.
### Also replaces the `FORMAT/GT` field of Y chromosome variants with "." for 
### all female samples (i.e., makes it a missing genotype call).

############################################################################


# accept VCF and manifest filepaths as command-line arguments
args <- commandArgs(trailingOnly = TRUE)
# check if VCF was provided
if ( length(args) < 2 ) {
  stop("Input VCF and manifest file are required",
       call. = FALSE)
}

vcfpath <- args[1]
manifestpath <- args[2]


############################################################################

### Replace full SV sequences in REF and ALT with <INS> or <DEL>

vcf <- fread(cmd = paste("zcat", vcfpath, "| grep -v '^#'"), 
             header = FALSE) %>%
  # if ref field has > 1 character, the SV is a deletion relative to the reference
  .[nchar(V4) > 1, V4 := "<DEL>"] %>%
  # if alt field has > 1 character, the SV is an insertion relative to the reference
  .[nchar(V5) > 1, V5 := "<INS>"]
  

############################################################################

### Replace Y chromosome genotype calls in female samples with "."

# replace the first character of a sample's FORMAT string,
# i.e. the genotype call for a Y chr variant, with "." (missing genotype)
set_gt_uncalled <- function(original_gt) {
  uncalled_gt <- sub("^.", # "^" - beginning of string; "." - any character
                     ".",
                     original_gt)
  return(uncalled_gt)
}

# read manifest file, select only female samples
manifest <- fread(manifestpath,
                  header = FALSE, col.names = c("path", "cov", "sex", "pop", "idx")) %>%
  .[ sex == "female" ]

# get index #s of female samples from manifest
female_idxs <- manifest$idx %>%
  # add 9 to index #s because the first 9 columns of the VCF are not sample columns
  + 9 %>%
  paste("V", ., sep = "") # column names for female samples

# change genotypes for female samples at all Y chr variants
vcf[ V1 == "chrY" ,
     (female_idxs) := pbmclapply(.SD,
                                 function(x) set_gt_uncalled(x),
                                 mc.cores = getOption("mc.cores", 24L)),
                      .SDcols = female_idxs]


############################################################################

### Output edited VCF to a file called `noseq.vcf`,
### for conversion back to a regular VCF by adding on a header

fwrite(vcf, "noseq.vcf",
       sep = "\t", row.names = FALSE, col.names = FALSE)