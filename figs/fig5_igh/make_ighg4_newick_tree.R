library(data.table)

### Modify Ohana's selection hypothesis C matrix for visualization of 
### ancestry component relationships at the IGHG4 locus as a Newick graph.
###
### Afterwards, run:
### ohana/bin/convert cov2nwk <c_matrix> <nwk.out>
### ohana/bin/convert nwk2svg <nwk.out> <svg.out>

############################################################################

### DATA

# path to "neutral" C matrix inferred from downsampled variants
c_neutral_path <- "chr21_pruned_50_C.matrix"
# path to C matrix of ancestry component where selection occurred
c_selection_path <- "k8/c_matrices/chr21_pruned_50_C_p2.matrix"

############################################################################


id <- "IGHG4" # name of SV
max_steps <- 100 # max steps used for selscan
step <- 32 # steps for SV of interest (from `selscan` output)
pop_num <- 2 # ancestry component #
k <- 8 # k for admixture

# SV's distance from neutral matrix, as fraction of the max distance
linear_interp <- step/max_steps

# neutral C matrix
c_neut <- fread(c_neutral_path,
                skip = 1)
# selection hypothesis C matrix
c_sel <- fread(c_selection_path,
              skip = 1)

new_c <- copy(c_neut)

# function to generate selection C matrix for that SV
make_cmat <- function(p, step) {
  a <- c(p-1)
  # max distance in selection matrix
  max_dist <- c_sel[p-1, ..a]
  # distance in neutral matrix
  min_dist <- c_neut[p-1, ..a]
  diff <- max_dist - min_dist
  sel_dist <- min_dist + diff * linear_interp
  
  # modify appropriate cell of C matrix
  new_c[p-1, p-1] = sel_dist
  
  ### log scaling (if necessary)
  # this is overly complicated because the distances are between 0 and 1
  # so most log values will be negative
  log_c <- log(new_c)
  # add this value to C matrix so that its min is 0
  zero <- abs(min(log(new_c)))
  # to avoid having branches with length 0, make the minimum branch length
  # equal to the minimum branch length of the non-logged matrix
  old_min <- min(new_c)
  # generate logged C matrix
  new_c <- log_c + zero + old_min
  
  # write header for Ohana C matrix
  write(paste(k-1, k-1, sep = " "),
        file = paste0(id, "_sel_log.matrix"))
  # write matrix to file
  fwrite(new_c,
         paste0(id, "_sel_log.matrix"),
         sep = " ", append = TRUE)
}

make_cmat(pop_num, step)