################################################################################
# Bayesian Re-analysis of Postoperative Delirium Prevention Interventions
#
# Based on: Meng et al. BJA 2024; 133(3): 565-583
#
# Bayesian hierarchical meta-analysis with posterior superiority probabilities
################################################################################

# ==============================================================================
# 0. Setup
# ==============================================================================

required_packages <- c(
  "tidyverse", "bayesmeta", "metafor",
  "ggplot2", "gridExtra", "scales", "RColorBrewer"
)

install_if_missing <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
}
invisible(lapply(required_packages, install_if_missing))
invisible(lapply(required_packages, library, character.only = TRUE))

set.seed(2024)

# ==============================================================================
# 1. Data Loading and Preparation
# ==============================================================================

df_main <- read.csv("pod_trials_main.csv", stringsAsFactors = FALSE)
df_rob  <- read.csv("pod_trials_rob.csv", stringsAsFactors = FALSE)

df <- df_main %>%
  left_join(df_rob %>% select(study_id, rob_overall), by = "study_id")

df_usual <- df %>% filter(control == "Usual care")

# Intervention categories
primary_interventions <- c(
  "Dexmedetomidine", "Non-pharmacological programmes",
  "Nerve or regional block", "EEG monitoring or light anaesthesia",
  "Insomnia treatment", "Antipsychotics", "Haemodynamics", "Steroids"
)

sensitivity_interventions <- c(
  "Electrical stimulation", "Alternative analgesia", "Respiration"
)

# Log OR and SE with continuity correction
calc_effect_size <- function(data) {
  data %>%
    mutate(
      pod_int_adj  = ifelse(pod_int == 0 | pod_int == n_int, pod_int + 0.5, pod_int),
      pod_ctrl_adj = ifelse(pod_ctrl == 0 | pod_ctrl == n_ctrl, pod_ctrl + 0.5, pod_ctrl),
      n_int_adj    = ifelse(pod_int == 0 | pod_int == n_int, n_int + 1, n_int),
      n_ctrl_adj   = ifelse(pod_ctrl == 0 | pod_ctrl == n_ctrl, n_ctrl + 1, n_ctrl),
      log_or = log((pod_int_adj / (n_int_adj - pod_int_adj)) /
                   (pod_ctrl_adj / (n_ctrl_adj - pod_ctrl_adj))),
      se_log_or = sqrt(1/pod_int_adj + 1/(n_int_adj - pod_int_adj) +
                       1/pod_ctrl_adj + 1/(n_ctrl_adj - pod_ctrl_adj))
    )
}

df_usual <- calc_effect_size(df_usual)

# ==============================================================================
# 2. Bayesian Meta-analysis Functions
# ==============================================================================

run_bayesian_ma <- function(data,
                            mu_prior_mean = 0,
                            mu_prior_sd = 10,
                            tau_prior_scale = 0.5) {
  if (nrow(data) < 3) return(NULL)

  tau_prior <- function(t) dt(t / tau_prior_scale, df = 1) / tau_prior_scale

  tryCatch(
    bayesmeta(
      y      = data$log_or,
      sigma  = data$se_log_or,
      labels = paste0(data$authors, " (", data$year, ")"),
      mu.prior.mean = mu_prior_mean,
      mu.prior.sd   = mu_prior_sd,
      tau.prior      = tau_prior
    ),
    error = function(e) { message("  Error: ", e$message); NULL }
  )
}

calc_posterior_probs <- function(result, thresholds = c(1, 0.9, 0.8, 0.7)) {
  if (is.null(result)) return(NULL)

  probs <- sapply(thresholds, function(th) result$pposterior(mu = log(th)))
  names(probs) <- paste0("P_OR_lt_", thresholds)

  summary_stats <- list(
    pooled_or       = exp(result$summary["mean", "mu"]),
    or_95cri_lower  = exp(result$summary["95% lower", "mu"]),
    or_95cri_upper  = exp(result$summary["95% upper", "mu"]),
    tau             = result$summary["mean", "tau"],
    tau_95cri_lower = result$summary["95% lower", "tau"],
    tau_95cri_upper = result$summary["95% upper", "tau"]
  )

  c(summary_stats, as.list(probs))
}

calc_superiority_matrix <- function(results_list, n_samples = 5000) {
  interventions <- names(results_list)
  n_int <- length(interventions)
  mat <- matrix(NA, n_int, n_int,
                dimnames = list(interventions, interventions))

  posterior_samples <- lapply(interventions, function(int) {
    if (!is.null(results_list[[int]])) {
      results_list[[int]]$rposterior(n = n_samples)[, "mu"]
    }
  })
  names(posterior_samples) <- interventions

  for (i in 1:n_int) {
    for (j in 1:n_int) {
      if (i != j &&
          !is.null(posterior_samples[[interventions[i]]]) &&
          !is.null(posterior_samples[[interventions[j]]])) {
        mat[i, j] <- mean(posterior_samples[[interventions[i]]] <
                          posterior_samples[[interventions[j]]])
      }
    }
  }
  mat
}

# ==============================================================================
# 3. Main Analysis
# ==============================================================================

results_primary <- list()
summary_list <- list()

for (int in primary_interventions) {
  int_data <- df_usual %>% filter(intervention == int)
  if (nrow(int_data) >= 3) {
    result <- run_bayesian_ma(int_data)
    if (!is.null(result)) {
      results_primary[[int]] <- result
      summary_list[[int]] <- c(
        intervention = int,
        n_trials  = nrow(int_data),
        n_patients = sum(int_data$n_int) + sum(int_data$n_ctrl),
        calc_posterior_probs(result)
      )
    }
  }
}

summary_df <- bind_rows(lapply(summary_list, as.data.frame.list)) %>%
  mutate(across(where(is.numeric), ~ round(., 3)))

cat("\n--- Primary Analysis Summary ---\n")
print(summary_df)

# Pairwise superiority matrix
superiority_mat <- calc_superiority_matrix(results_primary)

# ==============================================================================
# 4. Visualization
# ==============================================================================

# --- 4.1 Forest Plot ---
create_forest_plot <- function(summary_df) {
  plot_data <- summary_df %>%
    arrange(pooled_or) %>%
    mutate(
      intervention = factor(intervention, levels = intervention),
      label = sprintf("%.2f [%.2f, %.2f]", pooled_or, or_95cri_lower, or_95cri_upper)
    )

  ggplot(plot_data, aes(x = pooled_or, y = intervention)) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "gray40", linewidth = 0.5) +
    geom_errorbarh(aes(xmin = or_95cri_lower, xmax = or_95cri_upper),
                   height = 0.25, color = "gray40", linewidth = 0.4) +
    geom_point(aes(size = n_trials), color = "#2166AC", shape = 16) +
    geom_text(aes(x = 1.6, label = label), hjust = 0, size = 3.2, color = "gray20") +
    scale_x_continuous(
      trans = "log", breaks = c(0.3, 0.5, 0.7, 1.0, 1.5),
      limits = c(0.2, 2.0), expand = expansion(mult = c(0.02, 0.25))
    ) +
    scale_size_continuous(range = c(2.5, 7), name = "No. of trials") +
    labs(x = "Odds Ratio (95% CrI)", y = NULL,
         title = "Posterior Estimates of Intervention Effects on POD",
         subtitle = "Bayesian hierarchical meta-analysis vs usual care") +
    theme_minimal(base_size = 11) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.major.x = element_line(color = "gray90", linewidth = 0.3),
      legend.position = "bottom",
      plot.title = element_text(size = 12, face = "bold"),
      plot.subtitle = element_text(size = 9, color = "gray50")
    )
}

print(create_forest_plot(summary_df))

# --- 4.2 Superiority Heatmap ---
create_heatmap <- function(mat, summary_df) {
  mat_long <- as.data.frame(as.table(mat * 100))
  names(mat_long) <- c("Row", "Column", "Probability")

  order_int <- summary_df %>% arrange(pooled_or) %>% pull(intervention)
  mat_long$Row    <- factor(mat_long$Row, levels = order_int)
  mat_long$Column <- factor(mat_long$Column, levels = order_int)

  ggplot(mat_long, aes(x = Column, y = Row, fill = Probability)) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(aes(label = ifelse(is.na(Probability), "", sprintf("%.0f%%", Probability))),
              size = 3.2, color = "gray20") +
    scale_fill_gradient2(
      low = "#2166AC", mid = "white", high = "#B2182B",
      midpoint = 50, na.value = "gray95", limits = c(0, 100),
      name = "P(Row > Col)"
    ) +
    labs(x = NULL, y = NULL,
         title = "Posterior Probability of Superiority",
         subtitle = "P(Row intervention is more effective than Column intervention)") +
    theme_minimal(base_size = 11) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
      panel.grid = element_blank(),
      plot.title = element_text(size = 12, face = "bold"),
      plot.subtitle = element_text(size = 9, color = "gray50")
    )
}

print(create_heatmap(superiority_mat, summary_df))

# --- 4.3 Threshold Probability Plot ---
create_threshold_plot <- function(summary_df) {
  threshold_data <- summary_df %>%
    select(intervention, starts_with("P_OR_lt_")) %>%
    pivot_longer(cols = starts_with("P_OR_lt_"),
                 names_to = "threshold", values_to = "probability") %>%
    mutate(threshold = as.numeric(gsub("P_OR_lt_", "", threshold)),
           probability = probability * 100)

  ggplot(threshold_data, aes(x = threshold, y = probability,
                             color = intervention, group = intervention)) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 2.5) +
    geom_hline(yintercept = c(80, 95), linetype = "dashed", color = "gray50", alpha = 0.5) +
    annotate("text", x = 0.68, y = 82, label = "80%", size = 3, color = "gray50") +
    annotate("text", x = 0.68, y = 97, label = "95%", size = 3, color = "gray50") +
    scale_x_reverse(breaks = c(1, 0.9, 0.8, 0.7)) +
    scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20)) +
    scale_color_brewer(palette = "Set2", name = "Intervention") +
    labs(x = "OR Threshold", y = "P(OR < Threshold) %",
         title = "Probability of Effect Below Various OR Thresholds",
         subtitle = "Higher values indicate stronger evidence for risk reduction") +
    theme_minimal(base_size = 11) +
    theme(
      legend.position = "right", panel.grid.minor = element_blank(),
      plot.title = element_text(size = 12, face = "bold"),
      plot.subtitle = element_text(size = 9, color = "gray50")
    )
}

print(create_threshold_plot(summary_df))

# ==============================================================================
# 5. Sensitivity Analysis: Risk of Bias
# ==============================================================================

run_subgroup_sensitivity <- function(df_usual, interventions, filter_col, filter_vals, label) {
  df_filtered <- df_usual %>% filter(.data[[filter_col]] %in% filter_vals)

  results <- list()
  for (int in interventions) {
    int_data <- df_filtered %>% filter(intervention == int)
    if (nrow(int_data) >= 3) {
      result <- run_bayesian_ma(int_data)
      if (!is.null(result)) {
        results[[int]] <- c(
          intervention = int, filter = label,
          n_trials  = nrow(int_data),
          n_patients = sum(int_data$n_int) + sum(int_data$n_ctrl),
          calc_posterior_probs(result)
        )
      }
    }
  }

  if (length(results) > 0) {
    bind_rows(lapply(results, as.data.frame.list))
  } else {
    NULL
  }
}

# Risk of Bias
rob_combined <- bind_rows(
  run_subgroup_sensitivity(df_usual, primary_interventions, "rob_overall", "Low", "Low risk only"),
  run_subgroup_sensitivity(df_usual, primary_interventions, "rob_overall", c("Low", "Some concerns"), "Low + Some concerns"),
  summary_df %>% mutate(filter = "All trials") %>% select(intervention, filter, n_trials, n_patients, everything())
) %>% mutate(across(where(is.numeric), ~ round(., 3)))

cat("\n--- Sensitivity: Risk of Bias ---\n")
print(rob_combined)

# ==============================================================================
# 6. Sensitivity Analysis: Region
# ==============================================================================

region_combined <- bind_rows(
  run_subgroup_sensitivity(df_usual, primary_interventions, "region", "China", "China only"),
  run_subgroup_sensitivity(df_usual, primary_interventions, "region",
                           c("USA and Canada", "Europe Australia and New Zealand", "Middle East", "Other"),
                           "Non-China"),
  summary_df %>% mutate(filter = "All regions") %>% select(intervention, filter, n_trials, n_patients, everything())
) %>% mutate(across(where(is.numeric), ~ round(., 3)))

cat("\n--- Sensitivity: Region ---\n")
print(region_combined)

# ==============================================================================
# 7. Sensitivity Analysis: Prior Distributions
# ==============================================================================

prior_configs <- list(
  Weakly_informative     = list(mu_sd = 10, tau_scale = 0.5),
  Moderately_informative = list(mu_sd = 2,  tau_scale = 0.3),
  Skeptical              = list(mu_sd = 1,  tau_scale = 0.2)
)

prior_results <- list()
for (prior_name in names(prior_configs)) {
  config <- prior_configs[[prior_name]]
  summary_list_p <- list()
  for (int in primary_interventions) {
    int_data <- df_usual %>% filter(intervention == int)
    if (nrow(int_data) >= 3) {
      result <- run_bayesian_ma(int_data,
                                mu_prior_sd = config$mu_sd,
                                tau_prior_scale = config$tau_scale)
      if (!is.null(result)) {
        summary_list_p[[int]] <- c(
          intervention = int, prior = prior_name,
          n_trials  = nrow(int_data),
          n_patients = sum(int_data$n_int) + sum(int_data$n_ctrl),
          calc_posterior_probs(result)
        )
      }
    }
  }
  if (length(summary_list_p) > 0) {
    prior_results[[prior_name]] <- bind_rows(lapply(summary_list_p, as.data.frame.list))
  }
}

prior_combined <- bind_rows(prior_results) %>%
  mutate(across(where(is.numeric), ~ round(., 3)))

cat("\n--- Sensitivity: Prior Distributions ---\n")
print(prior_combined)

# ==============================================================================
# 8. Key Findings
# ==============================================================================

cat("\n--- Key Findings ---\n")

high_prob <- summary_df %>% filter(P_OR_lt_1 > 0.95) %>% arrange(desc(P_OR_lt_1))
cat("\nInterventions with P(OR < 1) > 95%:\n")
if (nrow(high_prob) > 0) {
  cat(sprintf("  %s: %.1f%%\n", high_prob$intervention, high_prob$P_OR_lt_1 * 100))
} else cat("  None\n")

clinical <- summary_df %>% filter(P_OR_lt_0.8 > 0.80) %>% arrange(desc(P_OR_lt_0.8))
cat("\nInterventions with P(OR < 0.8) > 80% (clinically meaningful):\n")
if (nrow(clinical) > 0) {
  cat(sprintf("  %s: %.1f%%\n", clinical$intervention, clinical$P_OR_lt_0.8 * 100))
} else cat("  None\n")
