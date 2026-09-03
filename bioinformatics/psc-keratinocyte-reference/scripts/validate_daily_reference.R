#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(ggplot2))

set.seed(350035)

root <- normalizePath(".")
gse_dir <- file.path(root, "output", "psc_kc_reference", "GSE120107")
daily_dir <- file.path(root, "output", "psc_kc_reference", "daily_reference_D0_D35")
out_dir <- file.path(daily_dir, "validation")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

expr <- read.delim(
  gzfile(file.path(gse_dir, "normalized_log2cpm.tsv.gz")),
  check.names = FALSE,
  row.names = 1
)
meta <- read.csv(file.path(gse_dir, "sample_qc.csv"), stringsAsFactors = FALSE)
stage_reference <- read.csv(file.path(daily_dir, "daily_stage_reference.csv"),
                            stringsAsFactors = FALSE)

anchor_stage <- c("D0", "D4", "D7", "D15", "D30", "Primary_KC")
anchor_day <- c(D0 = 0, D4 = 4, D7 = 7, D15 = 15, D30 = 30, Primary_KC = 35)
meta$target_day <- unname(anchor_day[meta$stage])

core_genes <- c(
  "POU5F1", "NANOG", "TFAP2A", "KRT18", "KRT19",
  "TP63", "KRT5", "KRT14", "ITGA6",
  "KRT1", "IVL", "SPRR1B", "FLG", "LOR", "ABCA12"
)
expanded_genes <- c(core_genes, "DSG1", "TGM1", "COL1A1", "PAX6")
panels <- list(core = core_genes, expanded = expanded_genes)

build_template <- function(train_samples, genes, stages = anchor_stage) {
  genes <- intersect(genes, rownames(expr))
  anchors <- sapply(stages, function(stage) {
    samples <- intersect(train_samples, meta$sample[meta$stage == stage])
    if (!length(samples)) stop("No training sample for stage: ", stage)
    rowMeans(expr[genes, samples, drop = FALSE])
  })
  colnames(anchors) <- stages
  scale_mean <- rowMeans(anchors)
  scale_sd <- apply(anchors, 1, sd)
  valid <- is.finite(scale_sd) & scale_sd > 0.25
  anchors <- anchors[valid, , drop = FALSE]
  scale_mean <- scale_mean[valid]
  scale_sd <- scale_sd[valid]
  daily <- t(apply(anchors, 1, function(values) {
    approx(unname(anchor_day[stages]), values, xout = 0:35,
           method = "linear", rule = 2)$y
  }))
  daily_z <- sweep(sweep(daily, 1, scale_mean, "-"), 1, scale_sd, "/")
  colnames(daily_z) <- paste0("D", 0:35)
  list(genes = rownames(anchors), mean = scale_mean, sd = scale_sd,
       daily_z = daily_z)
}

predict_one <- function(sample_name, train_samples, genes, n_boot = 500,
                        stages = anchor_stage) {
  template <- build_template(train_samples, genes, stages = stages)
  query <- expr[template$genes, sample_name]
  query_z <- (query - template$mean) / template$sd
  residual <- sweep(template$daily_z, 1, query_z, "-")
  rmse <- sqrt(colMeans(residual^2))
  correlations <- apply(template$daily_z, 2, function(x) {
    suppressWarnings(cor(x, query_z, method = "pearson"))
  })
  predicted_day <- which.min(rmse) - 1
  near <- which(rmse <= min(rmse) + 0.10) - 1

  boot_days <- replicate(n_boot, {
    idx <- sample(seq_along(template$genes), replace = TRUE)
    boot_rmse <- sqrt(colMeans(residual[idx, , drop = FALSE]^2))
    which.min(boot_rmse) - 1
  })
  ci <- unname(quantile(boot_days, c(0.025, 0.5, 0.975), type = 1))
  top <- order(rmse)[1:3] - 1
  stage_row <- stage_reference[stage_reference$day == predicted_day, ]

  data.frame(
    sample = sample_name,
    actual_stage = meta$stage[match(sample_name, meta$sample)],
    actual_day_or_state_anchor = meta$target_day[match(sample_name, meta$sample)],
    predicted_day = predicted_day,
    absolute_error_days = abs(predicted_day - meta$target_day[match(sample_name, meta$sample)]),
    bootstrap_median_day = ci[2],
    bootstrap_lower_day = ci[1],
    bootstrap_upper_day = ci[3],
    near_optimal_day_min = min(near),
    near_optimal_day_max = max(near),
    best_rmse = rmse[predicted_day + 1],
    best_correlation = correlations[predicted_day + 1],
    second_best_day = top[2],
    third_best_day = top[3],
    predicted_state = stage_row$expected_biological_state,
    predicted_confidence_class = stage_row$confidence,
    n_genes = length(template$genes),
    stringsAsFactors = FALSE
  )
}

prediction_rows <- list()
idx <- 1
all_samples <- meta$sample

for (panel_name in names(panels)) {
  genes <- panels[[panel_name]]

  for (query in all_samples) {
    result <- predict_one(query, setdiff(all_samples, query), genes)
    result$validation_method <- "leave_one_sample_out"
    result$panel <- panel_name
    result$training_replicate <- NA_integer_
    prediction_rows[[idx]] <- result
    idx <- idx + 1
  }

  for (train_rep in sort(unique(meta$replicate))) {
    train_samples <- meta$sample[meta$replicate == train_rep]
    test_samples <- meta$sample[meta$replicate != train_rep]
    for (query in test_samples) {
      result <- predict_one(query, train_samples, genes)
      result$validation_method <- "cross_replicate_trajectory"
      result$panel <- panel_name
      result$training_replicate <- train_rep
      prediction_rows[[idx]] <- result
      idx <- idx + 1
    }
  }

  # A harder interpolation test: remove every sample from one internal stage,
  # then ask the remaining flanking anchors to recover that missing stage.
  for (held_stage in c("D4", "D7", "D15", "D30")) {
    train_samples <- meta$sample[meta$stage != held_stage]
    test_samples <- meta$sample[meta$stage == held_stage]
    available_stages <- setdiff(anchor_stage, held_stage)
    for (query in test_samples) {
      result <- predict_one(query, train_samples, genes, stages = available_stages)
      result$validation_method <- "leave_one_stage_out"
      result$panel <- panel_name
      result$training_replicate <- NA_integer_
      prediction_rows[[idx]] <- result
      idx <- idx + 1
    }
  }
}

predictions <- do.call(rbind, prediction_rows)
predictions <- predictions[, c(
  "validation_method", "panel", "training_replicate", "sample", "actual_stage",
  "actual_day_or_state_anchor", "predicted_day", "absolute_error_days",
  "bootstrap_median_day", "bootstrap_lower_day", "bootstrap_upper_day",
  "near_optimal_day_min", "near_optimal_day_max", "best_rmse", "best_correlation",
  "second_best_day", "third_best_day", "predicted_state",
  "predicted_confidence_class", "n_genes"
)]
write.csv(predictions, file.path(out_dir, "daily_mapping_validation_predictions.csv"),
          row.names = FALSE)

metric_rows <- do.call(rbind, lapply(split(predictions,
                                           interaction(predictions$validation_method,
                                                       predictions$panel, drop = TRUE)),
                                    function(d) {
  data.frame(
    validation_method = unique(d$validation_method),
    panel = unique(d$panel),
    n_samples = nrow(d),
    mean_absolute_error_days = mean(d$absolute_error_days),
    median_absolute_error_days = median(d$absolute_error_days),
    within_2_days_fraction = mean(d$absolute_error_days <= 2),
    within_4_days_fraction = mean(d$absolute_error_days <= 4),
    bootstrap_interval_coverage = mean(
      d$actual_day_or_state_anchor >= d$bootstrap_lower_day &
        d$actual_day_or_state_anchor <= d$bootstrap_upper_day
    ),
    mean_best_correlation = mean(d$best_correlation),
    stringsAsFactors = FALSE
  )
}))
write.csv(metric_rows, file.path(out_dir, "daily_mapping_validation_metrics.csv"),
          row.names = FALSE)

plot_data <- subset(predictions, panel == "core")
method_labels <- c(
  cross_replicate_trajectory = "Cross-replicate trajectory",
  leave_one_sample_out = "Leave-one-sample-out",
  leave_one_stage_out = "Leave-one-stage-out interpolation"
)
plot_data$method_label <- unname(method_labels[plot_data$validation_method])
plot_data$method_label <- factor(
  plot_data$method_label,
  levels = unname(method_labels)
)
p <- ggplot(plot_data, aes(actual_day_or_state_anchor, predicted_day,
                           color = actual_stage)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
  geom_errorbar(aes(ymin = bootstrap_lower_day, ymax = bootstrap_upper_day),
                width = 0.5, alpha = 0.65) +
  geom_point(size = 3) +
  facet_wrap(~method_label, nrow = 1) +
  scale_x_continuous(breaks = unname(anchor_day)) +
  scale_y_continuous(breaks = seq(0, 35, 5), limits = c(0, 35)) +
  labs(
    title = "Daily reference self-validation",
    subtitle = "Core qPCR panel; D35 denotes a primary-keratinocyte state anchor",
    x = "Known public-data day/state anchor",
    y = "Predicted reference day",
    color = "Known stage"
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom")
ggsave(file.path(out_dir, "daily_mapping_validation.png"), p,
       width = 14, height = 6, dpi = 180)

capture.output(sessionInfo(), file = file.path(out_dir, "sessionInfo.txt"))
message("Validation written to: ", out_dir)
