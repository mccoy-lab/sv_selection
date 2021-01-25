library(EasyABC)
library(tidyverse)
library(cowplot)

runSLiM <- function(params) {
  seed <- as.integer(params[1])
  s_coef <- as.character(params[2])
  s_onset <- as.integer(round(params[3]))
  init_freq <- as.numeric(params[4])
  
  output <- system2("/home-net/home-4/rmccoy22@jhu.edu/code/slim/build/slim", c("-d", paste0("s_coef=", s_coef),
                                                                                "-d", paste0("s_onset=", s_onset),
                                                                                "-d", paste0("init_freq=", init_freq),
                                                                                "-s", seed, " /scratch/groups/rmccoy22/rmccoy22/sv_selection/slim/ighg4_growth.slim"), 
                    stdout = T, stderr = F)[16:20]
  freqs <- as.numeric(output)
  freqs[freqs == -Inf] <- 0
  freqs[freqs == Inf] <- 1
  if (init_freq < 0.5) freqs[is.na(freqs)] <- 0
  if (init_freq >= 0.5) freqs[is.na(freqs)] <- 1
  message(paste0("s = ", s_coef, ", ", 
                 "s_onset = ", s_onset, ", ", 
                 "init_freq = ", init_freq, ", ",
                 "AFS = ", freqs[1], ";", freqs[2], ";", 
                 freqs[3], ";", freqs[4], ";", freqs[5]))
  return(freqs)
}

# Set up and run our ABC
prior <- list(c("unif", -0.01, 0.2), c("unif", 1, 1685), c("unif", 0, 1))
observed <- c(0.09, 0.06, 0.27, 0.88, 0.65)

# ABC_SLiM <- ABC_rejection(model = runSLiM, prior = prior, nb_simul = 100000,
#                           tol = 0.01, summary_stat_target = observed,
#                           use_seed = TRUE, n_cluster = 44, seed_count = 1, progress_bar = TRUE)

ABC_SLiM <- ABC_sequential(method = "Lenormand", model = runSLiM, prior = prior, nb_simul = 10000,
                           summary_stat_target = observed,
                           use_seed = TRUE, n_cluster = 44, seed_count = 1, progress_bar = TRUE)

load("~/Downloads/ABC_SLiM.RData")

sum(ABC_SLiM$param[, 1] * ABC_SLiM$weights) 
sum(ABC_SLiM$param[, 2] * ABC_SLiM$weights) 
sum(ABC_SLiM$param[, 3] * ABC_SLiM$weights)
results <- data.table(s_coef = ABC_SLiM$param[, 1], 
                      s_onset = ABC_SLiM$param[, 2],
                      init_freq = ABC_SLiM$param[, 3],
                      weights = ABC_SLiM$weights,
                      CEU = ABC_SLiM$stats[, 1],
                      JPT = ABC_SLiM$stats[, 2],
                      CHB = ABC_SLiM$stats[, 3],
                      CDX = ABC_SLiM$stats[, 4],
                      KHV = ABC_SLiM$stats[, 5])

results[, s_onset := ((1686 - s_onset) * 29) / 1000]

a <- ggplot(data = results, aes(x = s_coef)) +
  #geom_histogram(data = data.table(prior = runif(1e7, min = as.numeric(prior[[1]][2]), max = as.numeric(prior[[1]][3]))),
  #               aes(x = prior, y = (..ncount..) * (ABC_SLiM$nsim / 50)), bins = 50) + 
  geom_histogram(bins = 60) +
  xlim(as.numeric(prior[[1]][2]), as.numeric(prior[[1]][3])) +
  ylab("Density") +
  xlab(expression(italic("s"))) + 
  theme_bw() +
  theme(panel.grid = element_blank())
  
b <- ggplot(data = results, aes(x = s_onset)) +
  #geom_histogram(data = data.table(prior = runif(ABC_SLiM$nsim, min = as.numeric(prior[[2]][2]), max = as.numeric(prior[[2]][3]))),
  #               aes(x = prior), bins = 50) + 
  geom_histogram(bins = 60) +
  #xlim(as.numeric(prior[[2]][2]), as.numeric(prior[[2]][3])) +
  xlim(0, 12) +
  ylab("Density") +
  xlab(expression(italic("T")["adaptive"]*" (KY)")) + 
  theme_bw() +
  theme(panel.grid = element_blank())

c <- ggplot(data = results, aes(x = init_freq)) +
  #geom_histogram(data = data.table(prior = runif(ABC_SLiM$nsim, min = as.numeric(prior[[3]][2]), max = as.numeric(prior[[3]][3]))),
  #               aes(x = prior), bins = 50) + 
  geom_histogram(bins = 60) +
  xlim(as.numeric(prior[[3]][2]), as.numeric(prior[[3]][3])) +
  ylab("Density") +
  xlab(expression(italic("p")[0])) + 
  theme_bw() +
  theme(panel.grid = element_blank())

d <- ggplot(results, aes(x = s_coef, y = s_onset)) +
  geom_point(size = 0.1) +
  xlim(as.numeric(prior[[1]][2]), as.numeric(prior[[1]][3])) +
  ylab(expression(italic("T")["adaptive"]*" (KY)")) + 
  xlab(expression(italic("s"))) + 
  theme_bw() +
  theme(panel.grid = element_blank())

cowplot::plot_grid(a, b, c, d, nrow = 2, 
                   align = "hv", axis = "l",
                   labels = c("B.", "C.", "D.", "E."))
