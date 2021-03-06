library(data.table)

### USAGE: Rscript /make_c_matrix.R <C.matrix file> <k>
### Generate covariance matrices to test for selection in one specific
### ancestry component with Ohana.
###
### <C.matrix file>: output from Ohana `nemeco`
### <k>: number of ancestry components, k, used to run `qpas`

############################################################################


# read in command line arguments
args <- commandArgs(trailingOnly = TRUE)
filepath <- args[1]
k <- args[2]
outpath <- paste0("k", k, "/c_matrices/")

# genome-wide "neutral" C matrix
cmat <- fread(filepath, skip = 1)

### C matrix for ancestry component 1
# write header for C matrix
write(paste(k, k, sep = " "),
      file = paste0(outpath, "chr21_pruned_50_C_p1.matrix"))
# write rest of C matrix for ancestry component 1:
# generated by adding a scalar (10) to all elements of the genome-wide matrix
# in accordance with Ohana's recommendation: https://github.com/jade-cheng/ohana/wiki/selscan
fwrite(cmat + 10,
       paste0(outpath, "chr21_pruned_50_C_p1.matrix"),
       sep = " ", append = TRUE)

### C matrix for ancestry components 2-k
# generated by adding a scalar (10) to the (k-1,k-1) element of the genome-wide matrix
# in accordance with Ohana's recommendation: https://github.com/jade-cheng/ohana/wiki/selscan
make_cmat <- function(p) {
  a <- c(p-1)
  val <- cmat[p-1, ..a]
  cmat[p-1, p-1] = val + 10
  
  # write header for C matrix
  write(paste(k, k, sep=" "),
        file = paste0(outpath, "chr21_pruned_50_C_p", p, ".matrix"))
  fwrite(cmat, 
         paste0(outpath, "chr21_pruned_50_C_p", p, ".matrix"),
         sep = " ", append = TRUE)
}
# apply to all ancestry components
invisible(lapply(2:k,
                 function(x) make_cmat(x)))
