# run_lme_metabolomics.R
# Linear mixed-effects models per metabolite
# Fixed effects: group (AX vs CV) + sex
# Random effect: parentage (random intercept)
# Also reports unadjusted group means and mean difference: mean CV - mean AX
# Inputs:
#   - merged_metabolomics.csv: CSV where the first column is metabolite ID/name and remaining columns are sample IDs
#   - SuppTableS5.xlsx: first column are sample IDs; must contain parentage, sex, and group columns
# Output:
#   - lme_results_all_group_plus_sex_fixed_with_mean_difference.csv in the working directory

suppressPackageStartupMessages({
  req <- c("readr","readxl","dplyr","tidyr","stringr","nlme","purrr")
  to_install <- setdiff(req, rownames(installed.packages()))
  if (length(to_install) > 0) install.packages(to_install, repos = "https://cloud.r-project.org")
  lapply(req, require, character.only = TRUE)
})

# ---- Config (edit paths if needed) ----
merged_path <- "merged_metabolomics.csv"
supp_path   <- "SuppTableS5.xlsx"
out_path    <- "lme_results_all_group_plus_sex_fixed_with_mean_difference.csv"

# ---- Helpers ----
find_group_col <- function(df) {
  name_hits <- which(tolower(names(df)) %in% c("group","treatment","condition"))
  if (length(name_hits) > 0) return(names(df)[name_hits[1]])

  for (nm in names(df)) {
    vals <- unique(na.omit(as.character(df[[nm]])))
    vals_upper <- toupper(vals)
    if (any(vals_upper %in% c("AX","CV"))) return(nm)
  }

  stop("Could not identify a 'group' column (AX/CV) in SuppTableS5.")
}

find_parentage_col <- function(df) {
  hit <- which(grepl("parentage", names(df), ignore.case = TRUE))
  if (length(hit) > 0) return(names(df)[hit[1]])
  stop("Could not find a 'parentage' column in SuppTableS5.")
}

find_sex_col <- function(df) {
  name_hits <- which(tolower(names(df)) %in% c("sex","gender"))
  if (length(name_hits) > 0) return(names(df)[name_hits[1]])

  hit <- which(grepl("sex", names(df), ignore.case = TRUE))
  if (length(hit) > 0) return(names(df)[hit[1]])

  stop("Could not find a 'sex' column in SuppTableS5.")
}

standardize_sex <- function(x) {
  x <- toupper(trimws(as.character(x)))
  x[x %in% c("M", "MALE", "MALES")] <- "M"
  x[x %in% c("F", "FEMALE", "FEMALES")] <- "F"
  x
}

read_supp_table <- function(path) {
  sheets <- tryCatch(readxl::excel_sheets(path), error = function(e) NULL)
  if (is.null(sheets) || length(sheets) == 0) {
    stop("No sheets found in SuppTableS5.xlsx.")
  }

  sheet_to_use <- if ("SuppTableS5" %in% sheets) "SuppTableS5" else sheets[1]
  supp <- readxl::read_excel(path, sheet = sheet_to_use)

  if (ncol(supp) < 4) {
    stop("SuppTableS5 must have sample_id plus metadata columns including parentage, sex, and group.")
  }

  sample_col <- names(supp)[1]
  parentage_col <- find_parentage_col(supp)
  group_col <- find_group_col(supp)
  sex_col <- find_sex_col(supp)

  supp2 <- supp |>
    dplyr::transmute(
      sample_id = as.character(.data[[sample_col]]),
      parentage = as.character(.data[[parentage_col]]),
      group = as.character(.data[[group_col]]),
      sex = as.character(.data[[sex_col]])
    ) |>
    dplyr::mutate(
      group = toupper(trimws(group)),
      sex = standardize_sex(sex),
      parentage = as.character(parentage)
    )

  supp2 <- dplyr::filter(
    supp2,
    group %in% c("AX","CV"),
    !is.na(sex),
    sex != ""
  )

  if (nrow(supp2) == 0) {
    stop("No rows in SuppTableS5 with group AX/CV and non-missing sex after processing.")
  }

  supp2
}

extract_coef <- function(tt, pattern) {
  rn <- rownames(tt)
  row <- rn[grepl(pattern, rn)]

  if (length(row) == 0) {
    return(list(
      estimate = NA_real_,
      std_error = NA_real_,
      df = NA_real_,
      t_value = NA_real_,
      p_value = NA_real_
    ))
  }

  row <- row[1]

  df_col <- intersect(colnames(tt), c("DF","denDF"))
  p_col <- intersect(colnames(tt), c("p-value","p-value:","p.value"))

  list(
    estimate = unname(tt[row, "Value"]),
    std_error = unname(tt[row, "Std.Error"]),
    df = if (length(df_col) == 0) NA_real_ else unname(tt[row, df_col[1]]),
    t_value = unname(tt[row, "t-value"]),
    p_value = if (length(p_col) == 0) NA_real_ else unname(tt[row, p_col[1]])
  )
}

compute_group_summaries <- function(df) {
  mean_AX <- mean(df$value[df$group == "AX"], na.rm = TRUE)
  mean_CV <- mean(df$value[df$group == "CV"], na.rm = TRUE)

  if (is.nan(mean_AX)) mean_AX <- NA_real_
  if (is.nan(mean_CV)) mean_CV <- NA_real_

  mean_difference <- if (!is.na(mean_AX) && !is.na(mean_CV)) {
    mean_CV - mean_AX
  } else {
    NA_real_
  }

  list(
    mean_AX = mean_AX,
    mean_CV = mean_CV,
    mean_difference_CV_minus_AX = mean_difference
  )
}

fit_one_metabolite <- function(df) {
  # df must contain columns: value, group (AX/CV), sex, parentage
  # AX is the reference group.
  # Sex reference is chosen alphabetically by factor(), usually F before M if values are F/M.
  group_summaries <- compute_group_summaries(df)

  df <- df |>
    dplyr::mutate(
      group = factor(group, levels = c("AX","CV")),
      sex = factor(sex),
      parentage = factor(parentage)
    ) |>
    dplyr::filter(!is.na(value), !is.na(group), !is.na(sex), !is.na(parentage))

  res <- list(
    group = list(estimate = NA_real_, std_error = NA_real_, df = NA_real_, t_value = NA_real_, p_value = NA_real_),
    sex = list(estimate = NA_real_, std_error = NA_real_, df = NA_real_, t_value = NA_real_, p_value = NA_real_),
    mean_AX = group_summaries$mean_AX,
    mean_CV = group_summaries$mean_CV,
    mean_difference_CV_minus_AX = group_summaries$mean_difference_CV_minus_AX,
    sigma = NA_real_,
    AIC = NA_real_,
    BIC = NA_real_,
    n = nrow(df),
    n_parentage = dplyr::n_distinct(df$parentage),
    n_sex = dplyr::n_distinct(df$sex),
    converged = FALSE,
    fit_warning = NA_character_
  )

  if (dplyr::n_distinct(df$group) < 2) {
    res$fit_warning <- "Insufficient variation in group"
    return(res)
  }

  if (dplyr::n_distinct(df$sex) < 2) {
    res$fit_warning <- "Insufficient variation in sex"
    return(res)
  }

  if (dplyr::n_distinct(df$parentage) < 1) {
    res$fit_warning <- "Insufficient variation in parentage"
    return(res)
  }

  fit <- tryCatch(
    nlme::lme(
      value ~ group + sex,
      random = ~ 1 | parentage,
      data = df,
      na.action = na.omit,
      control = nlme::lmeControl(msMaxIter = 100, msMaxEval = 200, returnObject = TRUE)
    ),
    error = function(e) e,
    warning = function(w) w
  )

  if (inherits(fit, "warning")) {
    res$fit_warning <- paste("warning:", conditionMessage(fit))
    return(res)
  }

  if (inherits(fit, "error")) {
    res$fit_warning <- paste("error:", conditionMessage(fit))
    return(res)
  }

  sm <- summary(fit)
  tt <- sm$tTable

  # Group coefficient is CV vs AX because AX is the reference level.
  res$group <- extract_coef(tt, "^group")

  # Sex coefficient is the non-reference sex level vs the reference sex level.
  # With F/M, this will usually be sexM = M vs F.
  res$sex <- extract_coef(tt, "^sex")

  if (is.na(res$group$estimate)) {
    res$fit_warning <- "No group coefficient found (unexpected design matrix)."
    return(res)
  }

  res$sigma <- sm$sigma
  res$AIC <- AIC(fit)
  res$BIC <- BIC(fit)
  res$converged <- TRUE
  res
}

message("Reading data...")
merged <- readr::read_csv(merged_path, show_col_types = FALSE)

if (ncol(merged) < 2) {
  stop("merged_metabolomics.csv must have at least 2 columns: metabolite plus at least one sample.")
}

metabolite_col <- names(merged)[1]
merged <- dplyr::rename(merged, metabolite = dplyr::all_of(metabolite_col))

long_df <- merged |>
  tidyr::pivot_longer(
    cols = -metabolite,
    names_to = "sample_id",
    values_to = "value"
  ) |>
  dplyr::mutate(sample_id = as.character(sample_id))

supp <- read_supp_table(supp_path)

dat <- long_df |>
  dplyr::inner_join(supp, by = "sample_id") |>
  dplyr::mutate(
    value = suppressWarnings(as.numeric(value))
  ) |>
  dplyr::filter(!is.na(value), !is.na(group), !is.na(sex), !is.na(parentage))

if (nrow(dat) == 0) {
  stop("No overlapping sample_id rows between merged_metabolomics.csv and SuppTableS5 with numeric values, AX/CV group, and sex.")
}

message("Fitting lme models per metabolite...")
results_list <- dat |>
  dplyr::group_by(metabolite) |>
  dplyr::group_split()

results <- purrr::map_dfr(results_list, function(dfm) {
  met <- unique(dfm$metabolite)
  fitres <- fit_one_metabolite(dfm)

  dplyr::tibble(
    metabolite = met,

    mean_AX = fitres$mean_AX,
    mean_CV = fitres$mean_CV,
    mean_difference_CV_minus_AX = fitres$mean_difference_CV_minus_AX,

    estimate_groupCV_vs_AX_adjusted_for_sex = fitres$group$estimate,
    std_error_group = fitres$group$std_error,
    df_group = fitres$group$df,
    t_value_group = fitres$group$t_value,
    p_value_group = fitres$group$p_value,

    estimate_sex_nonreference_vs_reference = fitres$sex$estimate,
    std_error_sex = fitres$sex$std_error,
    df_sex = fitres$sex$df,
    t_value_sex = fitres$sex$t_value,
    p_value_sex = fitres$sex$p_value,

    sigma = fitres$sigma,
    AIC = fitres$AIC,
    BIC = fitres$BIC,
    n_obs = fitres$n,
    n_parentage = fitres$n_parentage,
    n_sex = fitres$n_sex,
    converged = fitres$converged,
    note = fitres$fit_warning
  )
})

results <- results |>
  dplyr::arrange(dplyr::if_else(is.na(p_value_group), Inf, p_value_group))

readr::write_csv(results, out_path)
message(sprintf("Done. Wrote %s with %d rows.", out_path, nrow(results)))
