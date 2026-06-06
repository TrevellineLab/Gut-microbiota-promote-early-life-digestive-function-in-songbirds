### FIGURE 2
### 
###
### Model for each organ:
### log10_OrganMass ~ log10_TarsusLength + Treatment * Sex + (1 | ParentagePair)
###
### 
### 
###
###
###

suppressPackageStartupMessages({
  library(dplyr)
  library(lme4)
  library(lmerTest)
  library(tibble)
  library(ggplot2)
})

### LOAD DATA
df <- read.csv("DigestiveOrgans.csv", stringsAsFactors = FALSE)

### ADD SUM OF SMALL INTESTINE AND LARGE INTESTINE
df <- df %>%
  mutate(
    EmptySIplusLIMassMg =
      suppressWarnings(as.numeric(as.character(EmptySIMassMg))) +
      suppressWarnings(as.numeric(as.character(EmptyLIMassMg)))
  )

set.seed(32)

### ORGAN MASSES TO ANALYZE
organ_vars <- c(
  "LiverMassMg",
  "EmptySIMassMg",
  "EmptyLIMassMg",
  "EmptySIplusLIMassMg",
  "PancreasMassMg",
  "EmptyGizzardMassMg",
  "HeartMassMg"
)

### OPTIONAL: NICER FACET LABELS
organ_labels <- c(
  LiverMassMg = "Liver",
  EmptySIMassMg = "Small intestine",
  EmptyLIMassMg = "Large intestine",
  EmptySIplusLIMassMg = "Small + large intestine",
  PancreasMassMg = "Pancreas",
  EmptyGizzardMassMg = "Gizzard",
  HeartMassMg = "Heart"
)

### FACET LAYOUT
### ncol = 5 gives:
### Row 1: Liver | Small + large intestine | blank | blank | blank
### Row 2: Small intestine | Large intestine | Pancreas | Gizzard | Heart
facet_levels <- c(
  "Liver",
  "Small + large intestine",
  " ",
  "  ",
  "   ",
  "Small intestine",
  "Large intestine",
  "Pancreas",
  "Gizzard",
  "Heart"
)

### OUTPUT PREFIX
out_prefix <- "Figure2_organ_mass_LOG10_ANCOVA_lmer_Treatment_by_Sex_with_SIplusLI_no_clutch_color_no_lines"

### CHECK REQUIRED COLUMNS
required_cols <- c(
  "AnimalID",
  "Treatment",
  "Sex",
  "ParentagePair",
  "TarsusLength",
  organ_vars
)

missing_cols <- setdiff(required_cols, names(df))

if (length(missing_cols) > 0) {
  stop(
    "The following required columns are missing from DigestiveOrgans.csv: ",
    paste(missing_cols, collapse = ", ")
  )
}

### REMOVE OUTLIER NL113 IF PRESENT
df.filtered <- df %>%
  filter(AnimalID != "NL113") %>%
  mutate(
    Treatment = factor(Treatment),
    Sex = factor(Sex),
    ParentagePair = factor(ParentagePair)
  )

### STORAGE
type3_anova_all <- data.frame()
plot_data_all <- data.frame()

### HELPERS
as_numeric_safely <- function(x) {
  suppressWarnings(as.numeric(as.character(x)))
}

p_to_stars <- function(p) {
  if (length(p) == 0 || is.na(p)) {
    return("ns")
  }
  if (p < 0.001) {
    return("***")
  } else if (p < 0.01) {
    return("**")
  } else if (p < 0.05) {
    return("*")
  } else {
    return("ns")
  }
}

### HELPER: FISHER COMBINED STATISTIC
fisher_statistic_only <- function(p_values) {
  p_values <- suppressWarnings(as.numeric(p_values))
  p_values <- p_values[!is.na(p_values)]

  ### Avoid log(0) if any p-value is reported as exactly zero
  p_values[p_values <= 0] <- .Machine$double.xmin

  if (length(p_values) == 0) {
    return(data.frame(
      N_p_values = 0,
      Fisher_statistic = NA_real_
    ))
  }

  data.frame(
    N_p_values = length(p_values),
    Fisher_statistic = -2 * sum(log(p_values))
  )
}

### HELPER: GET SIX-ORGAN TREATMENT P-VALUES FOR CURRENT TREATMENT LABELS
get_treatment_pvalues_for_permutation_fisher <- function(input_df, fisher_organs, organ_labels) {

  pval_rows <- data.frame()

  for (organ in fisher_organs) {

    dat <- input_df %>%
      mutate(
        OrganMass_raw = as_numeric_safely(.data[[organ]]),
        TarsusLength_raw = as_numeric_safely(TarsusLength)
      ) %>%
      filter(
        !is.na(OrganMass_raw),
        !is.na(TarsusLength_raw),
        !is.na(Treatment),
        !is.na(Sex),
        !is.na(ParentagePair),
        OrganMass_raw > 0,
        TarsusLength_raw > 0
      ) %>%
      mutate(
        log10_OrganMass = log10(OrganMass_raw),
        log10_TarsusLength = log10(TarsusLength_raw),
        Treatment = droplevels(factor(Treatment)),
        Sex = droplevels(factor(Sex)),
        ParentagePair = droplevels(factor(ParentagePair)),
        Organ = organ,
        Organ_label = organ_labels[organ]
      )

    if (nrow(dat) < 4) next
    if (nlevels(dat$Treatment) < 2) next
    if (nlevels(dat$Sex) < 2) next
    if (nlevels(dat$ParentagePair) < 2) next

    full_fit <- tryCatch(
      lmerTest::lmer(
        log10_OrganMass ~ log10_TarsusLength + Treatment * Sex + (1 | ParentagePair),
        data = dat,
        REML = FALSE,
        control = lmerControl(
          optimizer = "bobyqa",
          optCtrl = list(maxfun = 100000)
        )
      ),
      error = function(e) e
    )

    if (inherits(full_fit, "error")) next

    anova_result <- tryCatch(
      {
        a <- as.data.frame(anova(full_fit, type = 3))
        a <- a %>% rownames_to_column(var = "Term")
        a$Organ <- organ
        a$Organ_label <- organ_labels[organ]
        a
      },
      error = function(e) e
    )

    if (inherits(anova_result, "error")) next

    p_col <- grep("^Pr\\(", names(anova_result), value = TRUE)
    if (length(p_col) == 0) next

    treatment_p <- anova_result[[p_col[1]]][anova_result$Term == "Treatment"]

    if (length(treatment_p) == 1 && !is.na(treatment_p)) {
      pval_rows <- bind_rows(
        pval_rows,
        data.frame(
          Organ = organ,
          Organ_label = organ_labels[organ],
          Treatment_p = treatment_p,
          stringsAsFactors = FALSE
        )
      )
    }
  }

  pval_rows
}

### HELPER: PERMUTE TREATMENT LABELS WITHIN PARENTAGEPAIR
permute_treatment_within_parentagepair <- function(input_df) {
  input_df %>%
    group_by(ParentagePair) %>%
    mutate(
      Treatment = sample(Treatment, size = n(), replace = FALSE)
    ) %>%
    ungroup() %>%
    mutate(Treatment = factor(Treatment))
}

### MAIN LOOP
for (organ in organ_vars) {

  cat("\n\n########################################\n")
  cat("ANALYSIS FOR:", organ, "\n")
  cat("########################################\n")

  dat <- df.filtered %>%
    mutate(
      OrganMass_raw = as_numeric_safely(.data[[organ]]),
      TarsusLength_raw = as_numeric_safely(TarsusLength)
    ) %>%
    filter(
      !is.na(OrganMass_raw),
      !is.na(TarsusLength_raw),
      !is.na(Treatment),
      !is.na(Sex),
      !is.na(ParentagePair),
      OrganMass_raw > 0,
      TarsusLength_raw > 0
    ) %>%
    mutate(
      log10_OrganMass = log10(OrganMass_raw),
      log10_TarsusLength = log10(TarsusLength_raw),
      Treatment = droplevels(factor(Treatment)),
      Sex = droplevels(factor(Sex)),
      ParentagePair = droplevels(factor(ParentagePair)),
      Organ = organ,
      Organ_label = organ_labels[organ],
      Sex_label = case_when(
        grepl("^m", tolower(as.character(Sex))) ~ "Male",
        grepl("^f", tolower(as.character(Sex))) ~ "Female",
        TRUE ~ as.character(Sex)
      )
    )

  if (nrow(dat) < 4) next
  if (nlevels(dat$Treatment) < 2) next
  if (nlevels(dat$Sex) < 2) next
  if (nlevels(dat$ParentagePair) < 2) next

  ### ANCOVA mixed model using log10 organ mass and log10 tarsus length
  full_fit <- tryCatch(
    lmerTest::lmer(
      log10_OrganMass ~ log10_TarsusLength + Treatment * Sex + (1 | ParentagePair),
      data = dat,
      REML = FALSE,
      control = lmerControl(
        optimizer = "bobyqa",
        optCtrl = list(maxfun = 100000)
      )
    ),
    error = function(e) e
  )

  if (inherits(full_fit, "error")) next

  ### Type III ANOVA
  anova_result <- tryCatch(
    {
      a <- as.data.frame(anova(full_fit, type = 3))
      a <- a %>% rownames_to_column(var = "Term")
      a$Organ <- organ
      a$Organ_label <- organ_labels[organ]
      a %>% select(Organ, Organ_label, everything())
    },
    error = function(e) e
  )

  if (inherits(anova_result, "error")) next

  type3_anova_all <- bind_rows(type3_anova_all, anova_result)

  ### Extract p-values for plot annotations
  p_col <- grep("^Pr\\(", names(anova_result), value = TRUE)
  if (length(p_col) == 0) next

  treatment_p <- anova_result[[p_col[1]]][anova_result$Term == "Treatment"]
  interaction_p <- anova_result[[p_col[1]]][anova_result$Term == "Treatment:Sex"]

  dat$Treatment_sig <- p_to_stars(treatment_p)
  dat$TreatmentSex_sig <- p_to_stars(interaction_p)

  plot_data_all <- bind_rows(plot_data_all, dat)
}

### CUSTOM FACET LAYOUT
if (nrow(plot_data_all) > 0) {
  plot_data_all$Organ_label <- factor(
    plot_data_all$Organ_label,
    levels = facet_levels
  )
}

### WRITE TYPE III ANOVA OUTPUT
write.csv(
  type3_anova_all,
  paste0(out_prefix, "_type3_anova.csv"),
  row.names = FALSE
)

### PERMUTATION-CALIBRATED FISHER TEST ACROSS DIGESTIVE ORGANS
### This tests for an aggregate Treatment effect across digestive organs.
### EmptySIplusLIMassMg and HeartMassMg are excluded from the aggregate test.
### Treatment labels are permuted within ParentagePair.
### The null distribution therefore preserves clutch structure and treatment counts within clutch.

fisher_organs <- c(
  "LiverMassMg",
  "EmptySIMassMg",
  "EmptyLIMassMg",
  "PancreasMassMg",
  "EmptyGizzardMassMg"
)

### Observed Treatment p-values from the five digestive organs
observed_treatment_pvalues <- get_treatment_pvalues_for_permutation_fisher(
  input_df = df.filtered,
  fisher_organs = fisher_organs,
  organ_labels = organ_labels
)

write.csv(
  observed_treatment_pvalues,
  paste0(out_prefix, "_observed_Treatment_pvalues_5_digestive_organs_for_permutation_Fisher.csv"),
  row.names = FALSE
)

observed_fisher <- fisher_statistic_only(observed_treatment_pvalues$Treatment_p) %>%
  mutate(Effect = "Treatment") %>%
  select(Effect, N_p_values, Fisher_statistic)

write.csv(
  observed_fisher,
  paste0(out_prefix, "_observed_Fisher_statistic_Treatment_5_digestive_organs.csv"),
  row.names = FALSE
)

### Permutation null distribution


n_perm <- 100

permutation_results <- data.frame()

for (perm_i in seq_len(n_perm)) {

  cat("Permutation", perm_i, "of", n_perm, "\n")

  df.perm <- permute_treatment_within_parentagepair(df.filtered)

  perm_pvalues <- get_treatment_pvalues_for_permutation_fisher(
    input_df = df.perm,
    fisher_organs = fisher_organs,
    organ_labels = organ_labels
  )

  perm_fisher <- fisher_statistic_only(perm_pvalues$Treatment_p)

  permutation_results <- bind_rows(
    permutation_results,
    data.frame(
      Permutation = perm_i,
      Effect = "Treatment",
      N_p_values = perm_fisher$N_p_values,
      Permuted_Fisher_statistic = perm_fisher$Fisher_statistic,
      stringsAsFactors = FALSE
    )
  )
}

observed_stat <- observed_fisher$Fisher_statistic[1]

permutation_results <- permutation_results %>%
  mutate(
    Observed_Fisher_statistic = observed_stat,
    Extreme_or_more_extreme = Permuted_Fisher_statistic >= Observed_Fisher_statistic
  )

empirical_p <- (sum(permutation_results$Extreme_or_more_extreme, na.rm = TRUE) + 1) /
  (sum(!is.na(permutation_results$Permuted_Fisher_statistic)) + 1)

permutation_summary <- data.frame(
  Effect = "Treatment",
  Permutation_scheme = "Treatment labels permuted within ParentagePair",
  N_permutations_requested = n_perm,
  N_permutations_successful = sum(!is.na(permutation_results$Permuted_Fisher_statistic)),
  N_observed_p_values = observed_fisher$N_p_values[1],
  Observed_Fisher_statistic = observed_stat,
  Empirical_p_value = empirical_p,
  stringsAsFactors = FALSE
)

write.csv(
  permutation_results,
  paste0(out_prefix, "_Permutation_null_Fisher_statistic_Treatment_5_digestive_organs_within_ParentagePair.csv"),
  row.names = FALSE
)

write.csv(
  permutation_summary,
  paste0(out_prefix, "_Permutation_summary_Fisher_statistic_Treatment_5_digestive_organs_within_ParentagePair.csv"),
  row.names = FALSE
)

### Plot observed Fisher statistic against permutation null distribution
permutation_plot_df <- permutation_results %>%
  filter(!is.na(Permuted_Fisher_statistic)) %>%
  mutate(
    Effect = factor(Effect, levels = c("Treatment"))
  )

observed_plot_df <- permutation_summary %>%
  mutate(
    Effect = factor(Effect, levels = c("Treatment")),
    Significance = case_when(
      Empirical_p_value < 0.001 ~ "***",
      Empirical_p_value < 0.01 ~ "**",
      Empirical_p_value < 0.05 ~ "*",
      TRUE ~ "ns"
    ),
    Significance_size = ifelse(Significance == "ns", 5, 7),
    Label = paste0(
      Significance,
      "\nPermutation p = ",
      signif(Empirical_p_value, 3)
    )
  )

p_perm <- ggplot(
    permutation_plot_df,
    aes(x = Effect, y = Permuted_Fisher_statistic)
  ) +
  geom_violin(
    width = 0.65,
    trim = TRUE,
    fill = "white",
    color = "black",
    linewidth = 0.6
  ) +
  geom_point(
    data = observed_plot_df,
    aes(x = Effect, y = Observed_Fisher_statistic),
    inherit.aes = FALSE,
    size = 3.2,
    color = "black"
  ) +
  geom_segment(
    data = observed_plot_df,
    aes(
      x = as.numeric(Effect) - 0.28,
      xend = as.numeric(Effect) + 0.28,
      y = Observed_Fisher_statistic,
      yend = Observed_Fisher_statistic
    ),
    inherit.aes = FALSE,
    linewidth = 1.1,
    color = "black"
  ) +
  geom_text(
    data = observed_plot_df,
    aes(
      x = Effect,
      y = Observed_Fisher_statistic,
      label = Label,
      size = Significance_size
    ),
    inherit.aes = FALSE,
    vjust = -0.8,
    fontface = "bold"
  ) +
  scale_size_identity() +
  scale_y_continuous(expand = expansion(mult = c(0.04, 0.14))) +
  labs(
    x = NULL,
    y = expression("Fisher statistic, " * X^2),
    title = "Observed aggregate Treatment statistic vs permutation null",
    subtitle = "Treatment labels permuted within ParentagePair; digestive organs only"
  ) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 0, hjust = 0.5)
  )

ggsave(
  filename = paste0(out_prefix, "_Permutation_null_Fisher_statistic_Treatment_5_digestive_organs_within_ParentagePair.pdf"),
  plot = p_perm,
  width = 4,
  height = 5
)

### WRITE PLOTTING DATA
write.csv(
  plot_data_all,
  paste0(out_prefix, "_raw_organ_mass_plot_data.csv"),
  row.names = FALSE
)

### MAKE PDF PLOT
if (nrow(plot_data_all) > 0) {

  annotation_df <- plot_data_all %>%
    group_by(Organ_label) %>%
    summarise(
      Treatment_sig = dplyr::first(Treatment_sig),
      TreatmentSex_sig = dplyr::first(TreatmentSex_sig),
      y_min = min(OrganMass_raw, na.rm = TRUE),
      y_max = max(OrganMass_raw, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      y_range = ifelse(y_max > y_min, y_max - y_min, abs(y_max)),
      y_range = ifelse(y_range == 0, 1, y_range),

      ### Force the comparison bar clearly above both violins/points
      bracket_y = y_max + 0.055 * y_range,
      bracket_tick_y = y_max + 0.025 * y_range,
      label_y = y_max + 0.095 * y_range,

      Treatment_sig_size = ifelse(Treatment_sig == "ns", 5, 8),
      TreatmentSex_sig_size = ifelse(TreatmentSex_sig == "ns", 5, 8)
    )

  p <- ggplot(
      plot_data_all,
      aes(x = Treatment, y = OrganMass_raw)
    ) +

    ### Narrow violin distributions
    geom_violin(
      width = 0.55,
      trim = TRUE,
      fill = "white",
      color = "black",
      linewidth = 0.5
    ) +

    ### Thick median bar
    stat_summary(
      fun = median,
      fun.min = median,
      fun.max = median,
      geom = "crossbar",
      width = 0.35,
      linewidth = 1.25,
      color = "black"
    ) +

    ### Individual data points at the forefront
    geom_point(
      aes(shape = Sex_label),
      position = position_jitter(width = 0.06, height = 0, seed = 1),
      size = 2.0,
      alpha = 0.95,
      color = "grey40",
      stroke = 0.4
    ) +

    facet_wrap(~ Organ_label, scales = "free_y", ncol = 5, drop = FALSE) +

    scale_y_continuous(expand = expansion(mult = c(0.04, 0.13))) +

    ### Treatment significance: bracket over the two Treatment violins
    geom_segment(
      data = annotation_df,
      aes(x = 1, xend = 2, y = bracket_y, yend = bracket_y),
      inherit.aes = FALSE,
      linewidth = 0.6,
      color = "black"
    ) +

    geom_segment(
      data = annotation_df,
      aes(x = 1, xend = 1, y = bracket_tick_y, yend = bracket_y),
      inherit.aes = FALSE,
      linewidth = 0.6,
      color = "black"
    ) +

    geom_segment(
      data = annotation_df,
      aes(x = 2, xend = 2, y = bracket_tick_y, yend = bracket_y),
      inherit.aes = FALSE,
      linewidth = 0.6,
      color = "black"
    ) +

    geom_text(
      data = annotation_df,
      aes(x = 1.5, y = label_y, label = Treatment_sig, size = Treatment_sig_size),
      inherit.aes = FALSE,
      vjust = 0.5,
      fontface = "bold"
    ) +

    scale_size_identity() +

    labs(
      x = "Treatment",
      y = "Raw organ mass, mg",
      shape = "Sex"
    ) +

    theme_bw() +
    theme(
      panel.grid = element_blank(),
      strip.background = element_rect(fill = "white"),
      legend.position = "right",
      strip.text = element_text(size = 9),
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.margin = margin(t = 8, r = 8, b = 8, l = 8)
    ) +
    coord_cartesian(clip = "off")

  ggsave(
    filename = paste0(out_prefix, "_raw_organ_mass_violinplot_custom_two_row_layout_TREATMENT_BRACKET_NO_TOPRIGHT_stats_from_model.pdf"),
    plot = p,
    width = 14,
    height = 8.5
  )
}

cat("\n\nDONE.\n")
cat("Wrote:\n")
cat(paste0(out_prefix, "_type3_anova.csv\n"))
cat(paste0(out_prefix, "_observed_Treatment_pvalues_5_digestive_organs_for_permutation_Fisher.csv\n"))
cat(paste0(out_prefix, "_observed_Fisher_statistic_Treatment_5_digestive_organs.csv\n"))
cat(paste0(out_prefix, "_Permutation_null_Fisher_statistic_Treatment_5_digestive_organs_within_ParentagePair.csv\n"))
cat(paste0(out_prefix, "_Permutation_summary_Fisher_statistic_Treatment_5_digestive_organs_within_ParentagePair.csv\n"))
cat(paste0(out_prefix, "_Permutation_null_Fisher_statistic_Treatment_5_digestive_organs_within_ParentagePair.pdf\n"))
cat(paste0(out_prefix, "_raw_organ_mass_plot_data.csv\n"))
cat(paste0(out_prefix, "_raw_organ_mass_violinplot_custom_two_row_layout_TREATMENT_BRACKET_NO_TOPRIGHT_stats_from_model.pdf\n"))
