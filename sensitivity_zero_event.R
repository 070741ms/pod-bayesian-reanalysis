################################################################################
# POD Bayesian Re-analysis: Zero-Event Sensitivity Analysis
#
# Compares three approaches to handling sparse-cell (zero-event) trials:
#   (A) Primary   : Standard 0.5 continuity correction (Haldane-Anscombe)
#   (B) Sweeting  : Reciprocal-of-other-arm correction (Sweeting et al. 2004)
#   (C) Exclusion : Sparse-cell trials removed
#
# Reference: Sweeting MJ, Sutton AJ, Lambert PC. What to add to nothing?
#   Use and avoidance of continuity corrections in meta-analysis of sparse
#   data. Stat Med. 2004;23(9):1351-1375. doi:10.1002/sim.1761
#
# Prior specification (identical to primary analysis):
#   mu  ~ Normal(0, 10^2)   [very weakly informative]
#   tau ~ Half-Cauchy(0, 0.5)
#
# Inputs  : pod_trials_main.csv, pod_trials_rob.csv  (working directory)
# Outputs : Sensitivity_ZeroEvent/Tables/  and  Sensitivity_ZeroEvent/Figures/
################################################################################

# ==============================================================================
# 0.  SETUP
# ==============================================================================

required_packages <- c("tidyverse", "bayesmeta", "scales", "patchwork")

install_if_missing <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE))
    install.packages(pkg, repos = "https://cloud.r-project.org")
}
invisible(lapply(required_packages, install_if_missing))
suppressPackageStartupMessages(
  invisible(lapply(required_packages, library, character.only = TRUE))
)

set.seed(2024)

# Output directories
OUT_DIR <- "Sensitivity_ZeroEvent"
TABLES  <- file.path(OUT_DIR, "Tables")
FIGURES <- file.path(OUT_DIR, "Figures")
for (d in c(OUT_DIR, TABLES, FIGURES))
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)

# Analysis constants (identical to primary analysis)
MU_PRIOR_MEAN   <- 0
MU_PRIOR_SD     <- 10
TAU_PRIOR_SCALE <- 0.5
MIN_TRIALS      <- 3

# Intervention categories (ordered by primary pooled OR, most protective first)
PRIMARY_INTERVENTIONS <- c(
  "Haemodynamics",
  "Dexmedetomidine",
  "Antipsychotics",
  "Nerve or regional block",
  "Steroids",
  "Non-pharmacological programmes",
  "Insomnia treatment",
  "EEG monitoring or light anaesthesia"
)

# ==============================================================================
# 1.  DATA LOADING
# ==============================================================================

if (!file.exists("pod_trials_main.csv") || !file.exists("pod_trials_rob.csv"))
  stop(
    "Required CSV files not found in working directory.\n",
    "Please set working directory to the folder containing:\n",
    "  pod_trials_main.csv\n  pod_trials_rob.csv",
    call. = FALSE
  )

df_main <- read.csv("pod_trials_main.csv", stringsAsFactors = FALSE)
df_rob  <- read.csv("pod_trials_rob.csv",  stringsAsFactors = FALSE)

df_raw <- df_main %>%
  left_join(df_rob %>% select(study_id, rob_overall), by = "study_id") %>%
  filter(control == "Usual care")

cat(sprintf("Trials loaded: %d usual-care-controlled\n", nrow(df_raw)))

# ==============================================================================
# 2.  SPARSE-CELL INVENTORY
#
#     "Sparse" = zero events OR complete events (all-event) in either arm.
#     These are the cells that require continuity correction.
# ==============================================================================

df_raw <- df_raw %>%
  mutate(
    sparse_int  = (pod_int  == 0 | pod_int  == n_int),
    sparse_ctrl = (pod_ctrl == 0 | pod_ctrl == n_ctrl),
    sparse_any  = sparse_int | sparse_ctrl,
    sparse_both = sparse_int & sparse_ctrl
  )

zero_inventory <- df_raw %>%
  filter(intervention %in% PRIMARY_INTERVENTIONS) %>%
  group_by(intervention) %>%
  summarise(
    n_trials   = n(),
    n_sparse   = sum(sparse_any),
    pct_sparse = round(100 * n_sparse / n_trials, 1),
    .groups    = "drop"
  ) %>%
  arrange(desc(n_sparse))

cat("\nSparse-cell counts by intervention:\n")
print(zero_inventory, n = Inf)

n_sparse_primary <- sum(
  df_raw$sparse_any & df_raw$intervention %in% PRIMARY_INTERVENTIONS
)
n_total_primary <- nrow(
  df_raw %>% filter(intervention %in% PRIMARY_INTERVENTIONS)
)
cat(sprintf(
  "\nSparse-cell trials (primary categories): %d / %d (%.1f%%)\n",
  n_sparse_primary, n_total_primary,
  100 * n_sparse_primary / n_total_primary
))

write.csv(zero_inventory,
          file.path(TABLES, "zero_inventory_by_intervention.csv"),
          row.names = FALSE)

# ==============================================================================
# 3.  EFFECT-SIZE CALCULATION — THREE METHODS
# ==============================================================================

# --------------------------------------------------------------------------
# Method A: Standard 0.5 continuity correction (Haldane-Anscombe)
#
#   For sparse trials, add 0.5 to events and 0.5 to non-events in both arms:
#     ri = pod_int  + 0.5,  ni = n_int  + 1
#     rc = pod_ctrl + 0.5,  nc = n_ctrl + 1
#   Non-sparse trials are unchanged.
# --------------------------------------------------------------------------
apply_cc05 <- function(data) {
  data %>%
    mutate(
      ri        = if_else(sparse_any, pod_int  + 0.5, as.double(pod_int)),
      ni        = if_else(sparse_any, n_int    + 1.0, as.double(n_int)),
      rc        = if_else(sparse_any, pod_ctrl + 0.5, as.double(pod_ctrl)),
      nc        = if_else(sparse_any, n_ctrl   + 1.0, as.double(n_ctrl)),
      log_or    = log((ri / (ni - ri)) / (rc / (nc - rc))),
      se_log_or = sqrt(1/ri + 1/(ni - ri) + 1/rc + 1/(nc - rc))
    )
}

# --------------------------------------------------------------------------
# Method B: Sweeting reciprocal-of-other-arm correction
#
#   Correction is proportional to relative arm size:
#     delta_I = n_I / (n_I + n_C)    [correction for intervention arm]
#     delta_C = n_C / (n_I + n_C)    [correction for control arm]
#
#   Applied as:
#     ri = pod_int  + delta_I,  ni = n_int  + 2 * delta_I
#     rc = pod_ctrl + delta_C,  nc = n_ctrl + 2 * delta_C
#
#   When arm sizes are equal, delta_I = delta_C = 0.5 (identical to Method A).
#   When arm sizes differ, the larger arm receives a smaller per-cell addition,
#   reducing bias from the fixed-0.5 rule.
#
#   Reference: Sweeting et al. Stat Med. 2004;23(9):1351-1375.
# --------------------------------------------------------------------------
apply_sweeting <- function(data) {
  data %>%
    mutate(
      delta_I   = n_int  / (n_int + n_ctrl),
      delta_C   = n_ctrl / (n_int + n_ctrl),
      ri        = if_else(sparse_any, pod_int  + delta_I,     as.double(pod_int)),
      ni        = if_else(sparse_any, n_int    + 2 * delta_I, as.double(n_int)),
      rc        = if_else(sparse_any, pod_ctrl + delta_C,     as.double(pod_ctrl)),
      nc        = if_else(sparse_any, n_ctrl   + 2 * delta_C, as.double(n_ctrl)),
      log_or    = log((ri / (ni - ri)) / (rc / (nc - rc))),
      se_log_or = sqrt(1/ri + 1/(ni - ri) + 1/rc + 1/(nc - rc))
    )
}

# --------------------------------------------------------------------------
# Method C: Zero-event exclusion
#
#   All sparse-cell trials are removed. No continuity correction is applied.
#   Raw cell counts are used directly for non-sparse trials.
# --------------------------------------------------------------------------
apply_exclusion <- function(data) {
  data %>%
    filter(!sparse_any) %>%
    mutate(
      ri        = as.double(pod_int),
      ni        = as.double(n_int),
      rc        = as.double(pod_ctrl),
      nc        = as.double(n_ctrl),
      log_or    = log((ri / (ni - ri)) / (rc / (nc - rc))),
      se_log_or = sqrt(1/ri + 1/(ni - ri) + 1/rc + 1/(nc - rc))
    )
}

# Method registry
METHODS <- list(
  A = list(
    label    = "0.5 correction (primary)",
    short    = "CC 0.5",
    fn       = apply_cc05,
    linetype = "solid",
    color    = "#2166AC",
    shape    = 16L
  ),
  B = list(
    label    = "Sweeting correction",
    short    = "Sweeting",
    fn       = apply_sweeting,
    linetype = "dashed",
    color    = "#E08214",
    shape    = 17L
  ),
  C = list(
    label    = "Zero-event exclusion",
    short    = "Exclusion",
    fn       = apply_exclusion,
    linetype = "dotted",
    color    = "#B2182B",
    shape    = 15L
  )
)

# ==============================================================================
# 4.  BAYESIAN HIERARCHICAL META-ANALYSIS
# ==============================================================================

# Core bayesmeta wrapper — identical prior specification to primary analysis
run_bayesmeta <- function(data) {
  if (nrow(data) < MIN_TRIALS) return(NULL)
  tau_prior <- function(t) dt(t / TAU_PRIOR_SCALE, df = 1) / TAU_PRIOR_SCALE
  tryCatch(
    bayesmeta(
      y             = data$log_or,
      sigma         = data$se_log_or,
      labels        = paste0(data$authors, " (", data$year, ")"),
      mu.prior.mean = MU_PRIOR_MEAN,
      mu.prior.sd   = MU_PRIOR_SD,
      tau.prior     = tau_prior
    ),
    error   = function(e) { warning(sprintf("bayesmeta error: %s", e$message)); NULL },
    warning = function(w) { message(sprintf("bayesmeta warning: %s", w$message)); NULL }
  )
}

# Extract posterior summaries
extract_posterior <- function(result, intervention, method_label,
                              n_trials, n_patients) {
  if (is.null(result)) return(NULL)
  thresh    <- c(1.0, 0.9, 0.8, 0.7)
  probs_pct <- sapply(thresh, function(th)
    result$pposterior(mu = log(th)) * 100
  )
  data.frame(
    intervention = intervention,
    method       = method_label,
    n_trials     = n_trials,
    n_patients   = n_patients,
    pooled_or    = exp(result$summary["mean",      "mu"]),
    or_lower     = exp(result$summary["95% lower", "mu"]),
    or_upper     = exp(result$summary["95% upper", "mu"]),
    tau_median   = result$summary["mean",      "tau"],
    tau_lower    = result$summary["95% lower", "tau"],
    tau_upper    = result$summary["95% upper", "tau"],
    P_OR_lt_1_0  = probs_pct[1],
    P_OR_lt_0_9  = probs_pct[2],
    P_OR_lt_0_8  = probs_pct[3],
    P_OR_lt_0_7  = probs_pct[4],
    stringsAsFactors = FALSE
  )
}

# Run all methods x all interventions
cat("\n=== Running Bayesian Meta-Analyses ===\n")

all_results <- list()

for (mkey in names(METHODS)) {
  m    <- METHODS[[mkey]]
  df_m <- m$fn(df_raw)

  cat(sprintf("\n  Method %s: %s  (n eligible = %d)\n",
              mkey, m$label, nrow(df_m)))

  for (int in PRIMARY_INTERVENTIONS) {
    int_data <- df_m %>% filter(intervention == int)

    if (nrow(int_data) < MIN_TRIALS) {
      cat(sprintf("    %-44s  SKIP (n=%d after filtering)\n",
                  int, nrow(int_data)))
      next
    }

    result <- run_bayesmeta(int_data)

    if (!is.null(result)) {
      row <- extract_posterior(
        result       = result,
        intervention = int,
        method_label = m$label,
        n_trials     = nrow(int_data),
        n_patients   = sum(int_data$n_int) + sum(int_data$n_ctrl)
      )
      all_results[[paste(mkey, int, sep = "__")]] <- row

      cat(sprintf(
        "    %-44s  n=%2d  OR=%5.3f [%5.3f, %5.3f]  P(OR<0.8)=%5.1f%%\n",
        int, nrow(int_data),
        row$pooled_or, row$or_lower, row$or_upper, row$P_OR_lt_0_8
      ))
    }
  }
}

# Combine into tidy data frame; preserve intervention order from PRIMARY_INTERVENTIONS
df_res <- bind_rows(all_results) %>%
  mutate(
    method = factor(method, levels = sapply(METHODS, `[[`, "label")),
    intervention = factor(intervention, levels = PRIMARY_INTERVENTIONS)
  ) %>%
  mutate(across(where(is.numeric), ~ round(., 4)))

# ==============================================================================
# 5.  TABLES
# ==============================================================================

cat("\n=== Writing Tables ===\n")

# Full results (all metrics)
write.csv(df_res,
          file.path(TABLES, "zero_sensitivity_full.csv"),
          row.names = FALSE)
cat("  zero_sensitivity_full.csv\n")

# Compact table — corresponds to Supplementary Table S7 Panel A
tbl_compact <- df_res %>%
  mutate(
    `OR (95% CrI)`  = sprintf("%.2f [%.2f, %.2f]", pooled_or, or_lower, or_upper),
    `tau (95% CrI)` = sprintf("%.2f [%.2f, %.2f]", tau_median, tau_lower, tau_upper)
  ) %>%
  select(
    Intervention   = intervention,
    Method         = method,
    Trials         = n_trials,
    `OR (95% CrI)`,
    `tau (95% CrI)`,
    `P(OR<1) %`    = P_OR_lt_1_0,
    `P(OR<0.8) %`  = P_OR_lt_0_8
  ) %>%
  arrange(Intervention, Method)

write.csv(tbl_compact,
          file.path(TABLES, "zero_sensitivity_compact.csv"),
          row.names = FALSE)
cat("  zero_sensitivity_compact.csv  [-> Supplementary Table S7 Panel A]\n")

cat("\n--- Compact table ---\n")
print(as.data.frame(tbl_compact))

# Deviation from primary analysis
df_primary_ref <- df_res %>%
  filter(method == METHODS$A$label) %>%
  select(intervention, or_ref = pooled_or, p08_ref = P_OR_lt_0_8)

tbl_delta <- df_res %>%
  filter(method != METHODS$A$label) %>%
  left_join(df_primary_ref, by = "intervention") %>%
  mutate(
    delta_OR  = round(pooled_or   - or_ref,    3),
    delta_P08 = round(P_OR_lt_0_8 - p08_ref,   2)
  ) %>%
  select(
    Intervention              = intervention,
    Method                    = method,
    Trials                    = n_trials,
    `Pooled OR`               = pooled_or,
    `delta OR vs primary`     = delta_OR,
    `P(OR<0.8) %`             = P_OR_lt_0_8,
    `delta P08 vs primary (pp)` = delta_P08
  ) %>%
  arrange(Intervention, Method)

write.csv(tbl_delta,
          file.path(TABLES, "zero_sensitivity_delta.csv"),
          row.names = FALSE)
cat("  zero_sensitivity_delta.csv\n")

# Summary of maximum deviations
delta_summary <- tbl_delta %>%
  group_by(Method) %>%
  summarise(
    max_abs_delta_OR  = round(max(abs(`delta OR vs primary`),           na.rm = TRUE), 3),
    max_abs_delta_P08 = round(max(abs(`delta P08 vs primary (pp)`),     na.rm = TRUE), 1),
    .groups = "drop"
  )

write.csv(delta_summary,
          file.path(TABLES, "zero_deviation_summary.csv"),
          row.names = FALSE)
cat("  zero_deviation_summary.csv\n")
cat("\n--- Maximum deviations from primary analysis ---\n")
print(delta_summary)

# ==============================================================================
# 6.  FIGURES
# ==============================================================================

cat("\n=== Creating Figures ===\n")

# Shared publication theme
theme_pub <- function(base_size = 11) {
  theme_minimal(base_size = base_size) %+replace% theme(
    text             = element_text(colour = "grey20"),
    plot.title       = element_text(size = rel(1.1), face = "bold",
                                    hjust = 0, margin = margin(b = 4)),
    plot.subtitle    = element_text(size = rel(0.84), colour = "grey42",
                                    hjust = 0, margin = margin(b = 10)),
    plot.caption     = element_text(size = rel(0.72), colour = "grey55",
                                    hjust = 1, margin = margin(t = 8),
                                    lineheight = 1.3),
    axis.title       = element_text(size = rel(0.9)),
    axis.text        = element_text(size = rel(0.82)),
    legend.title     = element_text(size = rel(0.88), face = "bold"),
    legend.text      = element_text(size = rel(0.82)),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(colour = "grey92", linewidth = 0.3),
    plot.margin      = margin(12, 16, 12, 12)
  )
}

method_colors <- setNames(sapply(METHODS, `[[`, "color"),
                          sapply(METHODS, `[[`, "label"))
method_shapes <- setNames(sapply(METHODS, `[[`, "shape"),
                          sapply(METHODS, `[[`, "label"))

# --------------------------------------------------------------------------
# Figure 1: Three-method comparison forest plot
# --------------------------------------------------------------------------
N_INT     <- nlevels(df_res$intervention)
y_offsets <- c(-0.26, 0, 0.26)

plot_forest <- df_res %>%
  mutate(
    y_base = as.integer(intervention),
    m_idx  = as.integer(method),
    y_pos  = y_base + y_offsets[m_idx]
  )

int_labels_forest <- df_res %>%
  filter(method == METHODS$A$label) %>%
  arrange(intervention) %>%
  transmute(lab = sprintf("%s  (n = %d)", as.character(intervention), n_trials)) %>%
  pull(lab)

x_lo <- max(0.15, min(plot_forest$or_lower, na.rm = TRUE) * 0.80)
x_hi <- min(3.00, max(plot_forest$or_upper, na.rm = TRUE) * 1.15)

fig1 <- ggplot(plot_forest,
               aes(x = pooled_or, y = y_pos, colour = method, shape = method)) +
  annotate("rect",
           xmin = -Inf, xmax = Inf,
           ymin = seq(0.5, N_INT - 0.5, 2),
           ymax = seq(1.5, N_INT + 0.5, 2),
           fill = "grey96", alpha = 1) +
  geom_vline(xintercept = 1, linetype = "dashed",
             colour = "grey35", linewidth = 0.45) +
  geom_errorbarh(aes(xmin = or_lower, xmax = or_upper),
                 height = 0.10, linewidth = 0.55) +
  geom_point(size = 2.5) +
  scale_x_continuous(
    trans  = "log",
    breaks = c(0.2, 0.3, 0.4, 0.5, 0.7, 1.0, 1.5, 2.0),
    labels = c("0.2", "0.3", "0.4", "0.5", "0.7", "1.0", "1.5", "2.0"),
    limits = c(x_lo, x_hi)
  ) +
  scale_y_continuous(
    breaks = seq_len(N_INT),
    labels = int_labels_forest,
    expand = expansion(add = 0.6)
  ) +
  scale_colour_manual(values = method_colors, name = "Zero-event method") +
  scale_shape_manual(values  = method_shapes, name = "Zero-event method") +
  labs(
    x        = "Odds Ratio (95% Credible Interval, log scale)",
    y        = NULL,
    title    = "Sensitivity Analysis: Zero-Event Handling Methods",
    subtitle = paste(
      "Primary (0.5 correction), Sweeting correction, and zero-event exclusion.",
      "Interventions ordered by primary pooled OR (most protective first)."
    ),
    caption  = paste(
      "CC 0.5 = Haldane\u2013Anscombe continuity correction (primary analysis).",
      "Sweeting = reciprocal-of-other-arm correction (Sweeting et al., Stat Med 2004;23:1351-1375).",
      "Exclusion = all sparse-cell trials removed.",
      "Trial count n shown for primary analysis; exclusion arm may differ.",
      sep = "\n"
    )
  ) +
  theme_pub() +
  theme(
    legend.position  = "bottom",
    legend.direction = "horizontal",
    axis.text.y      = element_text(size = 8.5, hjust = 1, lineheight = 1.1)
  )

tryCatch({
  ggsave(file.path(FIGURES, "fig1_zero_sensitivity_forest.png"),
         fig1, width = 250, height = 210, units = "mm", dpi = 300, bg = "white")
  cat("  fig1_zero_sensitivity_forest.png\n")
}, error = function(e) message(sprintf("  [!] Figure 1 failed: %s", e$message)))

# --------------------------------------------------------------------------
# Figure 2: Delta plot — departure from primary analysis
#   Corresponds to Supplementary Figure S5
#   Panel (a): delta pooled OR
#   Panel (b): delta P(OR < 0.8)
# --------------------------------------------------------------------------
delta_plot <- df_res %>%
  filter(method != METHODS$A$label) %>%
  left_join(
    df_res %>%
      filter(method == METHODS$A$label) %>%
      select(intervention, or_primary = pooled_or, p08_primary = P_OR_lt_0_8),
    by = "intervention"
  ) %>%
  mutate(
    delta_OR  = pooled_or   - or_primary,
    delta_P08 = P_OR_lt_0_8 - p08_primary,
    intervention = factor(intervention, levels = levels(df_res$intervention))
  )

x_range_OR  <- max(max(abs(delta_plot$delta_OR),  na.rm = TRUE) * 1.25, 0.05)
x_range_P08 <- max(max(abs(delta_plot$delta_P08), na.rm = TRUE) * 1.25, 1.0)

p_delta_OR <- ggplot(delta_plot,
                     aes(x = delta_OR, y = intervention,
                         colour = method, shape = method)) +
  geom_vline(xintercept = 0, linetype = "dashed",
             colour = "grey35", linewidth = 0.45) +
  geom_point(size = 3.2) +
  scale_x_continuous(
    limits = c(-x_range_OR, x_range_OR),
    labels = label_number(accuracy = 0.01, style_positive = "plus")
  ) +
  scale_colour_manual(values = method_colors[2:3], name = NULL) +
  scale_shape_manual(values  = method_shapes[2:3], name = NULL) +
  labs(x = "\u0394 Pooled OR vs primary", y = NULL,
       title = "(a)  \u0394 Pooled Odds Ratio") +
  theme_pub() +
  theme(legend.position = "bottom",
        axis.text.y = element_text(size = 9))

p_delta_P08 <- ggplot(delta_plot,
                      aes(x = delta_P08, y = intervention,
                          colour = method, shape = method)) +
  geom_vline(xintercept = 0, linetype = "dashed",
             colour = "grey35", linewidth = 0.45) +
  geom_point(size = 3.2) +
  scale_x_continuous(
    limits = c(-x_range_P08, x_range_P08),
    labels = label_number(accuracy = 0.1, suffix = " pp", style_positive = "plus")
  ) +
  scale_colour_manual(values = method_colors[2:3], name = NULL) +
  scale_shape_manual(values  = method_shapes[2:3], name = NULL) +
  labs(x = "\u0394 P(OR < 0.8) vs primary (percentage points)", y = NULL,
       title = "(b)  \u0394 P(OR < 0.8)") +
  theme_pub() +
  theme(legend.position = "bottom",
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())

fig2 <- (p_delta_OR | p_delta_P08) +
  plot_annotation(
    title    = "Impact of Zero-Event Handling on Posterior Estimates",
    subtitle = paste(
      "Horizontal displacement from zero indicates departure from the primary",
      "analysis. Interventions without sparse-cell trials show \u0394 = 0."
    ),
    caption  = paste(
      "Sweeting = reciprocal-of-other-arm correction.",
      "Exclusion = sparse-cell trials removed.",
      sep = "  |  "
    ),
    theme = theme_pub()
  ) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

tryCatch({
  ggsave(file.path(FIGURES, "fig2_zero_sensitivity_delta.png"),
         fig2, width = 240, height = 185, units = "mm", dpi = 300, bg = "white")
  cat("  fig2_zero_sensitivity_delta.png  [-> Supplementary Figure S5]\n")
}, error = function(e) message(sprintf("  [!] Figure 2 failed: %s", e$message)))

# ==============================================================================
# 7.  SUMMARY
# ==============================================================================

cat("\n================================================================\n")
cat("  ANALYSIS COMPLETE\n")
cat("================================================================\n\n")
cat(sprintf("Output folder: %s/\n\n", OUT_DIR))
cat("Tables/\n")
cat("  zero_inventory_by_intervention.csv\n")
cat("  zero_sensitivity_full.csv\n")
cat("  zero_sensitivity_compact.csv       <- Supplementary Table S7 Panel A\n")
cat("  zero_sensitivity_delta.csv\n")
cat("  zero_deviation_summary.csv\n\n")
cat("Figures/\n")
cat("  fig1_zero_sensitivity_forest.png\n")
cat("  fig2_zero_sensitivity_delta.png    <- Supplementary Figure S5\n\n")
cat("Session info:\n")
print(sessionInfo())
