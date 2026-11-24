# 99_proj_func.R -------------------------------------------------------------
# Project functions. Loaded in each .qmd with source("99_proj_func.R").


# Alpha diversity -------------------------------------------------------------
# All functions take a numeric vector of read counts for ONE sample (one value
# per ASV) and return a single number, so they can be used inside summarise():
#
#   dat |>
#     group_by(sample_id) |>
#     summarise(shannon = calc_shannon(reads))
#
# Zero counts are ignored, because absent ASVs do not contribute to diversity.


# Check that the input is a valid vector of read counts
check_counts <- function(counts) {
  if (!is.numeric(counts)) {
    stop("counts must be numeric")
  }
  if (any(is.na(counts))) {
    stop("counts must not contain NA")
  }
  if (any(counts < 0)) {
    stop("counts must not be negative")
  }
  if (sum(counts) == 0) {
    stop("counts must contain at least one read")
  }
}


# Observed richness: number of ASVs with at least one read
calc_observed <- function(counts) {
  check_counts(counts)
  sum(counts > 0)
}


# Chao1 estimated richness (bias-corrected):
#   S_obs + F1 * (F1 - 1) / (2 * (F2 + 1))
# F1 = number of singletons (ASVs seen once)
# F2 = number of doubletons (ASVs seen twice)
calc_chao1 <- function(counts) {
  check_counts(counts)
  s_obs <- sum(counts > 0)
  f1 <- sum(counts == 1)
  f2 <- sum(counts == 2)
  s_obs + f1 * (f1 - 1) / (2 * (f2 + 1))
}


# Shannon index: H' = -sum(p * ln(p)), p = relative abundance of each ASV
calc_shannon <- function(counts) {
  check_counts(counts)
  p <- counts[counts > 0] / sum(counts)
  -sum(p * log(p))
}


# Simpson index (Gini-Simpson): 1 - sum(p^2)
# Probability that two random reads belong to different ASVs
calc_simpson <- function(counts) {
  check_counts(counts)
  p <- counts[counts > 0] / sum(counts)
  1 - sum(p^2)
}


# Pielou evenness: H' / ln(S_obs)
# From 0 (one ASV dominates) to 1 (all ASVs equally abundant).
# Undefined (NA) when only one ASV is present.
calc_pielou <- function(counts) {
  check_counts(counts)
  s_obs <- sum(counts > 0)
  if (s_obs < 2) {
    return(NA_real_)
  }
  calc_shannon(counts) / log(s_obs)
}
