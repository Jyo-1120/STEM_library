#!/usr/bin/env Rscript

root <- normalizePath(".")
gse_dir <- file.path(root, "output", "psc_kc_reference", "GSE120107")
validation_dir <- file.path(
  root, "output", "psc_kc_reference", "daily_reference_D0_D35", "validation"
)
dir.create(validation_dir, recursive = TRUE, showWarnings = FALSE)

meta <- read.csv(file.path(gse_dir, "sample_qc.csv"), stringsAsFactors = FALSE)
pred <- read.csv(
  file.path(validation_dir, "daily_mapping_validation_predictions.csv"),
  stringsAsFactors = FALSE
)

anchor_day <- c(D0 = 0, D4 = 4, D7 = 7, D15 = 15, D30 = 30, Primary_KC = 35)
meta$target_day_or_state_anchor <- unname(anchor_day[meta$stage])

# Each fold trains on a complete time course from one replicate and validates
# on the other complete time course. No sample from the validation replicate is
# used to build that fold's reference.
split_rows <- list()
idx <- 1
for (fold in c("fold_1", "fold_2")) {
  train_rep <- if (fold == "fold_1") 1 else 2
  validation_rep <- if (fold == "fold_1") 2 else 1
  for (i in seq_len(nrow(meta))) {
    role <- if (meta$replicate[i] == train_rep) "training" else "validation"
    split_rows[[idx]] <- data.frame(
      fold = fold,
      training_replicate = train_rep,
      validation_replicate = validation_rep,
      role = role,
      sample = meta$sample[i],
      stage = meta$stage[i],
      target_day_or_state_anchor = meta$target_day_or_state_anchor[i],
      stringsAsFactors = FALSE
    )
    idx <- idx + 1
  }
}
split_manifest <- do.call(rbind, split_rows)
write.csv(
  split_manifest,
  file.path(validation_dir, "internal_train_validation_split.csv"),
  row.names = FALSE
)

cv <- subset(
  pred,
  validation_method == "cross_replicate_trajectory" & panel == "core"
)
cv$fold <- ifelse(cv$training_replicate == 1, "fold_1", "fold_2")
fold_metrics <- do.call(rbind, lapply(split(cv, cv$fold), function(d) {
  data.frame(
    fold = unique(d$fold),
    training_replicate = unique(d$training_replicate),
    validation_replicate = ifelse(unique(d$training_replicate) == 1, 2, 1),
    training_samples = 6,
    validation_samples = nrow(d),
    mean_absolute_error_days = mean(d$absolute_error_days),
    median_absolute_error_days = median(d$absolute_error_days),
    maximum_absolute_error_days = max(d$absolute_error_days),
    within_2_days_fraction = mean(d$absolute_error_days <= 2),
    mean_best_correlation = mean(d$best_correlation),
    bootstrap_interval_coverage = mean(
      d$actual_day_or_state_anchor >= d$bootstrap_lower_day &
        d$actual_day_or_state_anchor <= d$bootstrap_upper_day
    ),
    stringsAsFactors = FALSE
  )
}))
write.csv(
  fold_metrics,
  file.path(validation_dir, "internal_cross_replicate_fold_metrics.csv"),
  row.names = FALSE
)

external <- data.frame(
  dataset = c("GSE147206", "GSE287810", "GSE155816", "GSE98483", "GSE73305"),
  modality = c("scRNA-seq", "bulk RNA-seq", "scRNA-seq", "bulk RNA-seq", "bulk RNA-seq"),
  biological_material = c(
    "PSC-derived skin organoid epithelial cells",
    "iPSC-derived keratinocytes",
    "primary cultured keratinocytes",
    "primary keratinocytes",
    "primary keratinocytes"
  ),
  available_states = c(
    "D6, D29, D48",
    "D29-D32 endpoint; CnT30 versus KSFM",
    "passage 2 and passage 5",
    "calcium differentiation D0, D2, D4, D7",
    "undifferentiated and calcium D1-D5"
  ),
  validation_role = c(
    "cell-type specificity and D6/D29 state agreement",
    "late iKC endpoint agreement",
    "basal-to-differentiation marker direction",
    "primary-KC maturation direction",
    "terminal maturation direction"
  ),
  permitted_claim = c(
    "broad stage support; not exact calendar-day accuracy",
    "D30-like endpoint support; not exact day transfer",
    "maturation direction only",
    "maturation direction only",
    "maturation direction only"
  ),
  stringsAsFactors = FALSE
)
write.csv(
  external,
  file.path(validation_dir, "external_validation_manifest.csv"),
  row.names = FALSE
)

message("Train/validation manifests written to: ", validation_dir)
