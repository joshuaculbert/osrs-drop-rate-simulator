# OSRS Drop Rate Simulator
# This script simulates drop rates in Old School RuneScape (OSRS) to analyse:
# 1. Kills to first drop (geometric distribution)
# 2. Lucky streaks (drops in a 100-kill window)
# 3. Clustering of drops (inter-drop times)
# Outputs include statistical summaries and visualisations.
# Author: Joshua Culbert, March 2025

# 1. Setup and Parameters ------------------------------------------------------
library('ggplot2')

n_trials <- 10000          # Number of simulation trials 
n_kills <- 1000            # Number of kills per trial for first drop analysis
n_kills_clustering <- 5000 # Number of kills per trial for clustering analysis
window <- 100              # Window of kills for lucky streak analysis
drop_rates <- c(8, 16, 32, 64, 128, 256, 512, 1048, 3000, 5000) # 1/8 to 1/5000
probabilities <- 1 / drop_rates # Success probabilities
set.seed(123) # For reproducibility

# Create figures directory if not found
if (!dir.exists('figures')) {
  dir.create('figures')
}

# 2. Simulation Function -------------------------------------------------------
simulate_drops <- function(n_trials, n_kills, probs) {
  # Simulate binomial drops for given trials, kills, and probabilities
  # Args:
  #   n_trials: Number of trials
  #   n_kills: Number of kills per trial
  #   probs: Vector of success probabilities
  # Returns: 3D array of drops (trials x kills x probabilities)
  return(array(rbinom(n_trials * n_kills * length(probs),
                      size = 1,
                      prob = rep(probs, each = n_trials * n_kills)),
               dim = c(n_trials, n_kills, length(probs))))
}

# 3. Simulation of Kills to First Drop -----------------------------------------
# Simulate kills to first drop using geometric distribution
first_drops <- matrix(rgeom(n_trials * length(probabilities),
                            prob = rep(probabilities, each = n_trials)) + 1,
                            # Add 1 to account for 0-based indexing
                      nrow = n_trials)

# Convert to data frame for plotting
df_first_drops <- data.frame(
  kills = as.vector(first_drops), 
  prob = rep(paste('1/', drop_rates, sep = ''), each = n_trials)
)

# 4. Statistical Analysis ------------------------------------------------------
compute_stats <- function(data, probs, drop_rates) {
  # Compute simulated and theoretical statistics for geometric distribution
  # Args: 
  #   data: Matrix of kills to first drop (trials x drop rates)
  #   probs: Success probabilities
  #   drop_rates: Drop rates (1/probs)
  # Returns: Data frame with means, variances, IQRs, and outlier counts
  mean.sim <- apply(data, MARGIN = 2, FUN = mean)
  var.sim <- apply(data, MARGIN = 2, FUN = var)
  IQR.sim <- apply(first_drops, MARGIN = 2, FUN = IQR)
  outlier_thresh.sim <- sapply(1:length(probs), function(i) {
    quantile(data[, i], probs = 0.75) + 1.5 * IQR(data[, i])
  })
  outlier_count.sim <- sapply(1:length(probs), function(i) {
    sum(data[, i] > outlier_thresh.sim[i])
  })
  IQR.theoretical <- qgeom(p = 0.75, prob = probs) - 
                     qgeom(p = 0.25, prob = probs)
  outlier_thresh.theoretical <- qgeom(p = 0.75, prob = probs) +
                                1.5 * IQR.theoretical
  return(data.frame(
    prob = paste('1/', drop_rates, sep = ''),
    simulated_mean = round(mean.sim, 2),
    theoretical_mean = drop_rates,
    simulated_variance = round(var.sim, 2),
    theoretical_variance = (1 - probs) / probs^2,
    simulated_IQR = round(IQR.sim, 2),
    theoretical_IQR = IQR.theoretical,
    simulated_outlier_threshold = round(outlier_thresh.sim),
    theoretical_outlier_threshold = round(outlier_thresh.theoretical),
    simulated_outliers = outlier_count.sim,
    theoretical_outliers = round((1 - pgeom(q = outlier_thresh.theoretical, 
                                            prob = probs)) * nrow(data))
  ))
}

# Compute and output statistics
df_stats <- compute_stats(first_drops, probabilities, drop_rates)
cat('\nStatistical Summary of Kills to First Drop:\n')
print(df_stats)

# 5. CDF Visualisation ---------------------------------------------------------
# Theoretical geometric density with capped range
df_theoretical <- data.frame()
for (i in 1:length(probabilities)) {
  max_k <- min(max(first_drops[, i]), 10 / probabilities[i])  # Cap at 10 / p
  k <- 1:max_k
  dens <- dgeom(k - 1, prob = probabilities[i])  # Geometric density
  cum_dens <- pgeom(k - 1, prob = probabilities[i]) # Geometric CDF
  temp <- data.frame(
    kills = k, 
    density = dens, 
    cum_density = cum_dens,
    prob = paste('1/', drop_rates[i], sep = '')
  )
  df_theoretical <- rbind(df_theoretical, temp)
}

# Breaks list for log10 x-axis
breaks_list <- setNames(
  lapply(1 / probabilities, function(mean_kills) {
    max_k <- min(max(first_drops[, which(1 / probabilities == mean_kills)]),
                 10 * mean_kills)
    round(10^seq(0, log10(max_k), length.out = 5))
  }),
  paste('1/', drop_rates, sep = '')
)

# Plot empirical CDF vs. theoretical CDF
p_cdf <- ggplot() +
  stat_ecdf(
    data = df_first_drops,
    aes(x = kills, colour = 'Empirical'),
    linewidth = 1
  ) +
  geom_line(
    data = df_theoretical,
    aes(x = kills, y = cum_density, colour = 'Theoretical'),
    linewidth = 1
  ) +
  facet_wrap(~ prob, scales = 'free_x', ncol = 3) +
  scale_x_log10(
    breaks = function(x) {
      prob_label <- unique(
        df_theoretical$prob[df_theoretical$kills >= min(x) &
                            df_theoretical$kills <= max(x)]
      )
      if (length(prob_label) != 1 || !prob_label %in% names(breaks_list)) {
        return(c(1, 10, 100, 1000, 10000))
      }
      breaks_list[[prob_label]]
    }
  ) +
  scale_colour_manual(
    name = 'CDF Type',
    values = c('Empirical' = 'blue', 'Theoretical' = 'red')
  ) +
  coord_cartesian(xlim = c(1, NA)) +
  labs(
    title = 'Simulated vs. Theoretical Kills to First Drop',
    x = 'Number of Kills (Log10 Scale)',
    y = 'Cumulative Probability'
  ) +
  theme_minimal() +
  theme(
    strip.text = element_text(size = 10, colour = "white"), 
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, colour = "white"),
    axis.text.y  = element_text(colour = "white"),
    axis.title = element_text(colour = "white"),
    plot.title = element_text(colour = "white"),
    legend.text = element_text(colour = "white"),
    legend.title = element_text(colour = "white"),
    legend.position = 'bottom',
    plot.background = element_rect(fill = "#33333380", colour = NA),
    panel.background = element_rect(fill = "transparent", colour = NA)
  )

# Save and display the CDF plot
ggsave('figures/drop_rate_cdf_plot.png', p_cdf, 
       width = 12, height = 8, dpi = 300, bg = "#33333380")
print(p_cdf)

# 6. Lucky Streaks Analysis ----------------------------------------------------
# Simulate drops for streak analysis
drops <- simulate_drops(n_trials, n_kills, probabilities)

# Compute streak counts (drops in a 100-kill window)
streak_counts <- apply(drops[, 1:window, ], c(1, 3), sum)

# Calculate probabilities of 2+ drops
empirical_probs <- colMeans(streak_counts >= 2)
theoretical_probs <- 1 - ppois(q = 1, lambda = window * probabilities)

# Create summary table for streak probabilities
df_streak_probs <- data.frame(
  prob = paste('1/', drop_rates, sep = ''),
  empirical_prob_2plus = round(empirical_probs, 4),
  theoretical_prob_2plus = round(theoretical_probs, 4)
)
cat('\nProbability of 2+ Drops in a 100-Kill Window:\n')
print(df_streak_probs)

# Create data frame for plotting
df_streaks <- data.frame(
  streaks = as.vector(streak_counts),
  prob = rep(paste('1/', drop_rates, sep = ''), each = n_trials)
)

# Plot distribution of drops in a 100-kill window
p_streaks <- ggplot(df_streaks, aes(x = streaks)) +
  geom_histogram(
    aes(y = after_stat(density), fill = 'Empirical'),
    bins = 20,  # Increase bins for better resolution
    alpha = 0.5
  ) +
  geom_vline(
    aes(xintercept = 2, colour = 'Streak Threshold'),
    linetype = 'dashed',
    linewidth = 1
  ) +
  facet_wrap(~ prob, scales = 'free_y', ncol = 3) +
  scale_fill_manual(
    name = NULL,
    values = c('Empirical' = 'blue')
  ) +
  scale_colour_manual(
    name = NULL,
    values = c('Streak Threshold' = 'red')
  ) +
  guides(
    fill = guide_legend(title = 'Legend'),
    colour = guide_legend(title = NULL)
  ) +
  labs(
    title = 'Distribution of Drops in 100 Kills',
    x = 'Number of Drops',
    y = 'Density'
  ) +
  theme_minimal() +
  theme(
    strip.text = element_text(size = 10, colour = "white"), 
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, colour = "white"),
    axis.text.y  = element_text(colour = "white"),
    axis.title = element_text(colour = "white"),
    plot.title = element_text(colour = "white"),
    legend.text = element_text(colour = "white"),
    legend.title = element_text(colour = "white"),
    legend.position = 'bottom',
    plot.background = element_rect(fill = "#33333380", colour = NA),
    panel.background = element_rect(fill = "transparent", colour = NA)
  )

# Save and display the streaks plot
ggsave('figures/lucky_streaks_plot.png', p_streaks, 
       width = 12, height = 8, dpi = 300, bg = "#33333380")
print(p_streaks)

# 7. Clustering Analysis -------------------------------------------------------
# Simulate drop times for clustering analysis
drop_times <- simulate_drops(n_trials, n_kills_clustering, probabilities)

# Compute inter-drop times and KS tests
# Note: Small jitter (uniform, 0-1) is added to inter-drop times to break ties
# for the KS test, which assumes continuous data. Warnings are suppressed as
# they do not affect conclusions (see README for details).
ks_pvalues <- numeric(length(probabilities))
inter_times_list <- vector('list', length(probabilities))
for (i in 1:length(probabilities)) {
  # Get indices of drops across all trials
  drops <- which(drop_times[, , i] == 1, arr.ind = TRUE)
  if (nrow(drops) > 1) {
    # Sort by trial and then by kill index
    drops <- drops[order(drops[, 1], drops[, 2]), ]
    # Compute inter-drop times
    inter_times <- diff(drops[, 2])
    # Identify trial boundaries (where trial number changes)
    trial_changes <- which(diff(drops[, 1]) != 0)
    # Keep only inter-drop times within the same trial
    valid_indices <- setdiff(seq_along(inter_times), trial_changes)
    inter_times <- inter_times[valid_indices]
    # Apply jitter to break ties
    if (length(inter_times) > 0) {
      inter_times_list[[i]] <- inter_times + runif(length(inter_times), 0, 1)
      # KS test for exponential distribution, suppress warnings
      ks_result <- suppressWarnings(
        ks.test(inter_times_list[[i]], 'pexp', rate = probabilities[i])
      )
      ks_pvalues[i] <- ks_result$p.value
    } else {
      inter_times_list[[i]] <- numeric(0)
      ks_pvalues[i] <- NA
    }
  } else {
    inter_times_list[[i]] <- numeric(0)
    ks_pvalues[i] <- NA
  }
}

# Compute dispersion indices
dispersion_indices <- numeric(length(probabilities))
for (i in 1:length(probabilities)) {
  drop_counts <- colSums(matrix(drop_times[1, , i], nrow = window))
  dispersion_indices[i] <- var(drop_counts) / mean(drop_counts)
}

# Create data frame for inter-drop times plotting
df_inter_times <- data.frame(
  times = unlist(lapply(inter_times_list,
                        function(x) if (length(x) == 0) NA else x)),
  prob = rep(paste('1/', drop_rates, sep = ''),
             times = lengths(inter_times_list))
)
df_inter_times <- df_inter_times[!is.na(df_inter_times$prob), ]

# Plot inter-drop times
p_inter_times <- ggplot(df_inter_times, aes(x = times)) +
  geom_histogram(
    aes(y = after_stat(density), fill = 'Empirical'),
    bins = 30,
    alpha = 0.5
  ) +
  facet_wrap(~ prob, scales = 'free', ncol = 3)

# Add theoretical exponential curves for each drop rate
for (i in seq_along(drop_rates)) {
  max_x <- max(df_inter_times$times[df_inter_times$prob == 
                                      paste('1/', drop_rates[i], sep = '')], 
               na.rm = TRUE)
  p_inter_times <- p_inter_times +
    stat_function(
      data = data.frame(
        prob = paste('1/', drop_rates[i], sep = '')
      ),
      aes(colour = 'Theoretical'),  # Map colour for legend
      fun = dexp,
      args = list(rate = probabilities[i]),
      xlim = c(0, max_x),
      linewidth = 1,
      inherit.aes = FALSE
    )
}

p_inter_times <- p_inter_times +
  scale_fill_manual(
    name = 'Legend',
    values = c('Empirical' = 'blue', 'Theoretical' = 'red')
  ) +
  scale_colour_manual(
    name = 'Legend',
    values = c('Empirical' = 'blue', 'Theoretical' = 'red')
  ) +
  guides(
    fill = guide_legend(title = 'Legend'),
    colour = guide_legend(title = NULL, override.aes = list(fill = NA))
  ) +
  labs(
    title = 'Inter-Drop Times vs. Exponential Distribution',
    x = 'Kills Between Drops',
    y = 'Density'
  ) +
  theme_minimal() +
  theme(
    strip.text = element_text(size = 10, colour = "white"), 
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, colour = "white"),
    axis.text.y  = element_text(colour = "white"),
    axis.title = element_text(colour = "white"),
    plot.title = element_text(colour = "white"),
    legend.text = element_text(colour = "white"),
    legend.title = element_text(colour = "white"),
    legend.position = 'bottom',
    plot.background = element_rect(fill = "#33333380", colour = NA),
    panel.background = element_rect(fill = "transparent", colour = NA)
  )

# Save and display the inter-drop times plot
ggsave('figures/inter_drop_times_plot.png', p_inter_times, 
       width = 12, height = 8, dpi = 300, bg = "#33333380")
print(p_inter_times)

# Output clustering results
df_clustering <- data.frame(
  prob = paste('1/', drop_rates, sep = ''),
  ks_pvalue = sapply(ks_pvalues, function(x)
                                 if (is.na(x) || x < 0.01) "< 0.01" 
                                 else round(x, 2)),
  dispersion_index = round(dispersion_indices, 2)
)
cat('\nClustering Analysis Results:')
print(df_clustering)