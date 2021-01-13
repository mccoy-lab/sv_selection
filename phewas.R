library(data.table)
library(tidyverse)
library(pbmcapply)
library(qqman)

setwd("~/Downloads/phewas/")

manifest_url <- "https://pan-ukb-us-east-1.s3.amazonaws.com/sumstats_release/phenotype_manifest.tsv.bgz"
download.file(manifest_url, dest = basename(manifest_url))
manifest <- fread(paste0("gzcat ", basename(manifest_url)))

eas_phewas <- manifest[grepl("EAS", pops)]

get_phewas <- function(eas_phewas, row_index, download_index = TRUE) {
  result = tryCatch({
    url <- eas_phewas[row_index,]$aws_link
    options(timeout = 3600)
    if (download_index == TRUE) {
      tabix_url <- eas_phewas[row_index,]$aws_link_tabix
      download.file(tabix_url, destfile = basename(tabix_url))
    }
    header <- c("chr",	"pos",	"ref",	"alt",	"af_meta",	"beta_meta",	"se_meta",	"pval_meta",	
                "pval_heterogeneity",	"af_AFR",	"af_AMR",	"af_CSA",	"af_EAS",	"af_EUR",	"af_MID",	
                "beta_AFR",	"beta_AMR",	"beta_CSA",	"beta_EAS",	"beta_EUR",	"beta_MID",	"se_AFR",	
                "se_AMR",	"se_CSA",	"se_EAS",	"se_EUR",	"se_MID",	"pval_AFR",	"pval_AMR",	"pval_CSA",	
                "pval_EAS",	"pval_EUR",	"pval_MID",	"low_confidence_AFR",	"low_confidence_AMR",	
                "low_confidence_CSA",	"low_confidence_EAS",	"low_confidence_EUR",	"low_confidence_MID")
    tabix_cmd <- paste("tabix", url, "14:106226552-106226552")
    phewas_line <- fread(cmd = tabix_cmd)
    if (nrow(phewas_line) == 0) {
      phewas_line <- setNames(data.table(matrix(nrow = 0, ncol = length(header))), header) %>%
        .[, description := eas_phewas[row_index,]$description]
      return(phewas_line)
    } else {
      setnames(phewas_line, header)
      phewas_line[, description := eas_phewas[row_index,]$description]
      return(phewas_line)
    }
  }, error = function(error_condition) {
    phewas_line <- setNames(data.table(matrix(nrow = 0, ncol = length(header))), header) %>%
      .[, description := eas_phewas[row_index,]$description]
    return(phewas_line)
  })
}

phewas_results <- rbindlist(pbmclapply(1:nrow(eas_phewas), 
                                       function(x) suppressWarnings(get_phewas(eas_phewas, x, download_index = FALSE)), 
                                       mc.cores = getOption("mc.cores", 6L)))

setorder(phewas_results, pval_EAS)
qq(phewas_results$pval_EAS)
