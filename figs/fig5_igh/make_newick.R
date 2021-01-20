library(data.table)

### Modify Ohana's selection hypothesis C matrices for visualization as Newick graphs.
###
### I assume that the "step" result output from selscan reflects the selected SV's
### fraction of the linear distance between the neutral and selection matrices.
###
### Afterwards, run:
### ohana/bin/convert cov2nwk <c_matrix> <nwk.out>
### ohana/bin/convert nwk2svg <nwk.out> <svg.out>
###
### currently this DOES NOT work for the first ancestry component because
### the selection matrix is different and I didn't bother implementing it


# setwd("/Users/syan/Documents/mccoy-lab/sv_selection/ohana")
setwd("/scratch/groups/rmccoy22/syan11/sv_selection/ohana")

id <- "IGHG4" # name of SV
max_steps <- 100 # max steps used for selscan
step <- 32 # steps for SV of interest
pop_num <- 1 # ancestry component #
k <- 8 # k for admixture

# SV's distance from neutral matrix, as fraction of the max distance
linear_interp <- step/max_steps

# neutral C matrix
c_neut <- fread("k8/chr21_pruned_50_C.matrix",
                skip = 1)
# selection hypothesis C matrix
c_sel <- fread(paste0("k8/c_matrices/chr21_pruned_50_C_p",
                     pop_num, ".matrix"),
              skip = 1)

new_c <- copy(c_neut)

# function to generate selection C matrix for that SV
make_cmat <- function(p, step) {
  a <- c(p)
  # max distance in selection matrix
  max_dist <- c_sel[p, ..a]
  # distance in neutral matrix
  min_dist <- c_neut[p, ..a]
  diff <- max_dist - min_dist
  sel_dist <- min_dist + diff * linear_interp
  
  # modify appropriate cell of C matrix
  new_c[p, p] = sel_dist
  
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
  write(paste(k-1, k-1, sep=" "),
        file = paste0(id, "_sel_log.matrix"))
  # write matrix to file
  fwrite(new_c,
         paste0(id, "_sel_log.matrix"),
         sep = " ", append = TRUE)
}

make_cmat(pop_num, step)