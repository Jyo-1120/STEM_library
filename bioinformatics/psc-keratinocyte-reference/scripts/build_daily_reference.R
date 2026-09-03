#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(ggplot2))

root <- normalizePath(".")
input_dir <- file.path(root, "output", "psc_kc_reference", "GSE120107")
out_dir <- file.path(root, "output", "psc_kc_reference", "daily_reference_D0_D35")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

stage_expr <- read.delim(
  file.path(input_dir, "stage_mean_log2cpm.tsv"),
  check.names = FALSE,
  row.names = 1
)
module_sample <- read.csv(file.path(input_dir, "module_scores.csv"))

anchor_stage <- c("D0", "D4", "D7", "D15", "D30", "Primary_KC")
anchor_day <- c(0, 4, 7, 15, 30, 35)
anchor_kind <- c(
  "observed_bulk", "observed_bulk", "observed_bulk",
  "observed_bulk", "observed_bulk", "primary_KC_state_anchor"
)

genes <- c(
  "POU5F1", "NANOG", "TFAP2A", "KRT18", "KRT19",
  "TP63", "KRT5", "KRT14", "ITGA6", "KRT1", "KRT10",
  "IVL", "DSG1", "TGM1", "SPRR1B", "FLG", "LOR", "ABCA12",
  "COL1A1", "PAX6"
)
genes <- intersect(genes, rownames(stage_expr))

days <- 0:35

protocol_phase <- function(day) {
  if (day == 0) return("PSC / differentiation start")
  if (day <= 4) return("RA+BMP4 ectoderm induction")
  if (day == 5) return("RA+BMP4 passage")
  if (day == 6) return("UCM:N2 transition + RA/BMP4")
  if (day == 7) return("N2 + RA/BMP4")
  if (day == 8) return("RA/BMP4 withdrawal; EGF start")
  if (day <= 13) return("EGF keratinocyte progenitor induction")
  if (day == 14) return("keratinocyte progenitor checkpoint / passage")
  if (day <= 22) return("EGF progenitor expansion")
  if (day == 23) return("CaCl2 maturation start")
  if (day <= 29) return("CaCl2 keratinocyte maturation")
  if (day == 30) return("mature iKC checkpoint")
  "late iKC maintenance / primary-KC-state interpolation"
}

biological_state <- function(day) {
  if (day == 0) return("pluripotent PSC")
  if (day <= 3) return("ectoderm induction")
  if (day <= 8) return("surface ectoderm / epithelial induction")
  if (day <= 13) return("early keratinocyte progenitor")
  if (day <= 17) return("keratinocyte progenitor")
  if (day <= 22) return("basal-like immature keratinocyte")
  if (day <= 29) return("calcium-responsive maturing keratinocyte")
  if (day == 30) return("immature/maturing iKC endpoint")
  "late-maturing keratinocyte-like state"
}

evidence_class <- function(day) {
  if (day %in% c(0, 4, 7, 15, 30)) return("observed bulk anchor")
  if (day %in% c(6, 29)) return("interpolated; cross-protocol scRNA support")
  if (day == 35) return("primary KC state anchor; not an observed Day 35")
  if (day > 30) return("late endpoint interpolation")
  if (day %in% c(5, 8, 14, 23)) return("protocol event; expression interpolated")
  "piecewise-linear interpolation"
}

confidence <- function(day) {
  if (day %in% c(0, 4, 7, 15, 30)) return("high")
  if (day %in% c(6, 29)) return("medium-high")
  if (day > 30) return("low")
  if (day %in% c(5, 8, 14, 23)) return("medium")
  "medium-low"
}

stage_daily <- data.frame(
  day = days,
  day_label = paste0("D", days),
  protocol_phase = vapply(days, protocol_phase, character(1)),
  expected_biological_state = vapply(days, biological_state, character(1)),
  evidence_class = vapply(days, evidence_class, character(1)),
  confidence = vapply(days, confidence, character(1)),
  stringsAsFactors = FALSE
)

stage_daily$protocol_event <- ""
stage_daily$protocol_event[stage_daily$day == 0] <- "Start RA + BMP4 in UCM"
stage_daily$protocol_event[stage_daily$day == 5] <- "Passage 1:2; ROCK inhibitor 24 h"
stage_daily$protocol_event[stage_daily$day == 6] <- "UCM:N2 = 1:1"
stage_daily$protocol_event[stage_daily$day == 7] <- "Switch completely to N2"
stage_daily$protocol_event[stage_daily$day == 8] <- "Withdraw RA/BMP4; add EGF"
stage_daily$protocol_event[stage_daily$day == 14] <- "Progenitor qPCR/IF checkpoint and passage"
stage_daily$protocol_event[stage_daily$day == 15] <- "Remove ROCK inhibitor"
stage_daily$protocol_event[stage_daily$day == 23] <- "Add 1.2 mM CaCl2"
stage_daily$protocol_event[stage_daily$day == 30] <- "Mature-marker checkpoint"
stage_daily$protocol_event[stage_daily$day == 35] <- "Maintenance; model endpoint only"

write.csv(stage_daily, file.path(out_dir, "daily_stage_reference.csv"), row.names = FALSE)

gene_rows <- lapply(genes, function(gene) {
  anchor_values <- as.numeric(stage_expr[gene, anchor_stage])
  predicted <- approx(anchor_day, anchor_values, xout = days, method = "linear", rule = 2)$y
  observed_range <- max(anchor_values) - min(anchor_values)
  relative <- if (observed_range > 0) {
    (predicted - min(anchor_values)) / observed_range
  } else {
    rep(0.5, length(predicted))
  }
  level <- cut(
    relative,
    breaks = c(-Inf, 1/3, 2/3, Inf),
    labels = c("low", "intermediate", "high"),
    right = FALSE
  )
  z <- if (sd(predicted) > 0) as.numeric(scale(predicted)) else rep(0, length(predicted))
  data.frame(
    day = days,
    day_label = paste0("D", days),
    gene = gene,
    expected_log2cpm = predicted,
    within_gene_z = z,
    relative_to_gene_range = relative,
    expected_level = as.character(level),
    is_observed_bulk_anchor = days %in% c(0, 4, 7, 15, 30),
    is_primary_KC_state_anchor = days == 35,
    stringsAsFactors = FALSE
  )
})
gene_daily <- do.call(rbind, gene_rows)
write.csv(gene_daily, file.path(out_dir, "daily_gene_reference_long.csv"), row.names = FALSE)

gene_wide <- reshape(
  gene_daily[, c("day", "gene", "within_gene_z")],
  idvar = "day", timevar = "gene", direction = "wide"
)
names(gene_wide) <- sub("within_gene_z\\.", "z_", names(gene_wide))
daily_wide <- merge(stage_daily, gene_wide, by = "day", all.x = TRUE, sort = TRUE)
write.csv(daily_wide, file.path(out_dir, "daily_reference_for_mapping.csv"), row.names = FALSE)

module_mean <- aggregate(score ~ stage + stage_order + module, module_sample, mean)
module_rows <- lapply(unique(module_mean$module), function(module_name) {
  d <- module_mean[module_mean$module == module_name, ]
  d <- d[match(anchor_stage, d$stage), ]
  predicted <- approx(anchor_day, d$score, xout = days, method = "linear", rule = 2)$y
  data.frame(
    day = days,
    day_label = paste0("D", days),
    module = module_name,
    expected_module_score = predicted,
    is_observed_bulk_anchor = days %in% c(0, 4, 7, 15, 30),
    is_primary_KC_state_anchor = days == 35,
    stringsAsFactors = FALSE
  )
})
module_daily <- do.call(rbind, module_rows)
write.csv(module_daily, file.path(out_dir, "daily_module_reference_long.csv"), row.names = FALSE)

phase_lines <- data.frame(
  x = c(5, 8, 14, 23, 30),
  label = c("Passage", "EGF", "Progenitor", "CaCl2", "D30 checkpoint")
)

p_module <- ggplot(module_daily, aes(day, expected_module_score, color = module)) +
  geom_line(linewidth = 1) +
  geom_point(data = subset(module_daily, is_observed_bulk_anchor), size = 2) +
  geom_vline(data = phase_lines, aes(xintercept = x), inherit.aes = FALSE,
             linetype = "dashed", color = "grey70") +
  scale_x_continuous(breaks = seq(0, 35, 5)) +
  labs(
    title = "D0-D35 daily module reference",
    subtitle = "Points: observed GSE120107 bulk anchors; D35: primary-KC state anchor",
    x = "Protocol day", y = "Expected module score", color = "Module"
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom")
ggsave(file.path(out_dir, "daily_module_reference.png"), p_module,
       width = 11, height = 7, dpi = 180)

focus_genes <- c("POU5F1", "TFAP2A", "KRT18", "TP63", "KRT5", "KRT14",
                 "KRT1", "IVL", "SPRR1B", "FLG", "LOR")
p_gene <- ggplot(subset(gene_daily, gene %in% focus_genes),
                 aes(day, within_gene_z, color = gene)) +
  geom_line(linewidth = 0.9) +
  geom_point(data = subset(gene_daily, gene %in% focus_genes & is_observed_bulk_anchor),
             size = 1.8) +
  geom_vline(data = phase_lines, aes(xintercept = x), inherit.aes = FALSE,
             linetype = "dashed", color = "grey75") +
  facet_wrap(~gene, scales = "free_y", ncol = 4) +
  scale_x_continuous(breaks = c(0, 7, 15, 23, 30, 35)) +
  labs(
    title = "Selected-marker daily reference",
    subtitle = "Piecewise-linear interpolation; D31-D35 has low calendar-day confidence",
    x = "Protocol day", y = "Within-gene standardized expression"
  ) +
  theme_bw(base_size = 10) +
  theme(legend.position = "none")
ggsave(file.path(out_dir, "daily_marker_reference.png"), p_gene,
       width = 11, height = 8, dpi = 180)

anchor_table <- data.frame(
  model_day = anchor_day,
  source_stage = anchor_stage,
  anchor_kind = anchor_kind,
  interpretation = c(
    "PSC baseline",
    "early ectoderm induction",
    "surface ectoderm / epithelial induction",
    "keratinocyte progenitor transition",
    "immature/maturing iKC endpoint",
    "primary keratinocyte-like state; not a measured protocol day"
  )
)
write.csv(anchor_table, file.path(out_dir, "model_anchors.csv"), row.names = FALSE)

message("Daily reference written to: ", out_dir)
