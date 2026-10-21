#!/usr/bin/env Rscript

# Recompose the v6.12 main figures from already-existing result tables.
# This is a plotting-only script: it does not fit models or recalculate statistics.

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
})

# --- Path configuration: EDIT these placeholder roots for your machine ----------
DATA_ROOT    <- "<DATA_ROOT>"      # raw GEO / supplementary files (see README "Inputs")
WORK_ROOT    <- "<WORK_ROOT>"       # author analysis workspace (intermediate result tables)
ARCHIVE_ROOT <- "<ARCHIVE_ROOT>"    # author manuscript output archive
OUTPUT_ROOT  <- "./output"          # results are written here (relative to repository root)

root <- ARCHIVE_ROOT
real_out_dir <- file.path(root, "outputs/manuscript_v6_12_optionB_2026-10-13")
out_dir <- file.path(OUTPUT_ROOT, "fig_genetics/script2_manuscript_v6_12")
fig_dir <- file.path(out_dir, "figures")
src_out <- file.path(out_dir, "figure_source_data")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(src_out, recursive = TRUE, showWarnings = FALSE)

v611 <- file.path(root, "outputs/manuscript_v6_11_definition_alignment_2026-10-11")
tab <- file.path(v611, "tables")
native <- file.path(v611, "source_data")

ink <- "#25343D"
blue <- "#28789A"
copper <- "#B86D4B"
slate <- "#687984"
midgray <- "#A9B5BC"
pale_blue <- "#EAF2F5"
pale_copper <- "#F7EEE9"
lightgray <- "#EEF1F3"

theme_paper <- function(base_size = 10) {
  theme_classic(base_size = base_size, base_family = "sans") +
    theme(
      text = element_text(colour = ink),
      axis.title = element_text(size = base_size, colour = ink),
      axis.text = element_text(size = base_size - 1, colour = ink),
      axis.line = element_line(linewidth = 0.38, colour = "#82919A"),
      axis.ticks = element_line(linewidth = 0.35, colour = "#82919A"),
      plot.title = element_text(size = base_size + 1.4, face = "bold", colour = ink,
                                margin = margin(b = 5)),
      plot.subtitle = element_text(size = base_size - 0.3, colour = slate,
                                   margin = margin(b = 7)),
      plot.caption = element_text(size = base_size - 1.2, colour = slate, hjust = 0),
      panel.grid.major.x = element_line(linewidth = 0.28, colour = "#E6EAED"),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "top",
      legend.title = element_blank(),
      legend.text = element_text(size = base_size - 0.7),
      legend.key.width = grid::unit(0.9, "lines"),
      plot.margin = margin(8, 12, 8, 8)
    )
}

read_tsv <- function(path) {
  read.delim(path, check.names = FALSE, stringsAsFactors = FALSE,
             na.strings = c("NA", "", "NaN"))
}

save_figure <- function(plot, name, width, height) {
  ggsave(file.path(fig_dir, paste0(name, ".png")), plot,
         width = width, height = height, units = "in", dpi = 600,
         device = ragg::agg_png, bg = "white", limitsize = FALSE)
  ggsave(file.path(fig_dir, paste0(name, ".pdf")), plot,
         width = width, height = height, units = "in",
         device = grDevices::pdf, bg = "white", useDingbats = FALSE,
         limitsize = FALSE)
}

q_label <- function(x) {
  ifelse(is.na(x), "NA", formatC(x, format = "g", digits = 2))
}

# Preserve exact input tables alongside the figures so each display is traceable.
copy_source <- function(source, name) {
  target <- file.path(src_out, name)
  if (!file.copy(source, target, overwrite = TRUE)) stop("Could not copy source: ", source)
  target
}

# Figure 1: uses the exact cohort table extracted from the v6.12 manuscript.
cohorts <- read_tsv(file.path(real_out_dir, "figure_source_data", "Figure1_design_source.tsv"))
get_cohort <- function(id) {
  x <- cohorts[cohorts$Dataset == id, , drop = FALSE]
  if (nrow(x) != 1) stop("Expected one Table 1 row for ", id)
  x
}
get_group_counts <- function(row) {
  z <- strsplit(as.character(row[["MDD / control"]]), "/", fixed = TRUE)[[1]]
  if (length(z) != 2) stop("Unexpected case/control counts: ", row[["MDD / control"]])
  paste0(trimws(z[1]), " MDD / ", trimws(z[2]), " control")
}
one_blood <- function(id) {
  x <- get_cohort(id)
  specimen <- as.character(x[["Specimen and platform"]])
  specimen <- gsub("Blood leukocyte RNA sequencing", "Leukocytes, RNA-seq", specimen, fixed = TRUE)
  counts <- strsplit(x[["MDD / control"]], "/", fixed = TRUE)[[1]]
  n_total <- sum(as.numeric(counts))
  paste0(id, "  |  ", specimen, "\n", get_group_counts(x), ", n=", n_total)
}
g987 <- one_blood("GSE98793")
g396 <- one_blood("GSE39653")
g197 <- one_blood("GSE19738")
g251 <- one_blood("GSE251778")
paired <- get_cohort("GSE251780")
brain_ids <- c("GSE53987", "GSE54563", "GSE80655")
stopifnot(all(brain_ids %in% cohorts$Dataset))

cards <- data.frame(
  xmin = c(0.35, 4.05, 7.75), xmax = c(3.55, 7.25, 10.95),
  ymin = c(2.45, 2.45, 2.45), ymax = c(6.65, 6.65, 6.65),
  fill = c(pale_copper, pale_blue, "#F1F3F4")
)
fig1 <- ggplot() +
  geom_rect(data = cards, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = cards$fill, colour = "#D6DEE2", linewidth = 0.55) +
  annotate("text", x = 1.95, y = 6.18, label = "1  EXPLORATORY SCREEN",
           fontface = "bold", size = 3.3, colour = ink) +
  annotate("text", x = 1.95, y = 5.26, label = g987,
           size = 2.75, lineheight = 1.17, colour = ink) +
  annotate("text", x = 1.95, y = 4.36, label = g396,
           size = 2.75, lineheight = 1.17, colour = ink) +
  annotate("text", x = 1.95, y = 3.17,
           label = "Broad GO-BP screening\nOXPHOS carried forward",
           fontface = "bold", size = 2.55, lineheight = 1.2, colour = copper) +
  annotate("text", x = 5.65, y = 6.18, label = "2  FOCUSED FOLLOW-UP",
           fontface = "bold", size = 3.3, colour = ink) +
  annotate("text", x = 5.65, y = 5.26, label = g197,
           size = 2.75, lineheight = 1.17, colour = ink) +
  annotate("text", x = 5.65, y = 4.36, label = g251,
           size = 2.75, lineheight = 1.17, colour = ink) +
  annotate("text", x = 5.65, y = 3.17,
           label = "Correlation-aware CAMERA and FRY\nParticipant-level scores",
           fontface = "bold", size = 2.55, lineheight = 1.2, colour = blue) +
  annotate("text", x = 9.35, y = 6.18, label = "3  MOLECULAR CONTEXT",
           fontface = "bold", size = 3.3, colour = ink) +
  annotate("text", x = 9.35, y = 5.25,
           label = paste0("Paired blood methylation\nGSE251780: ", get_group_counts(paired),
                          "\n147 participants; subset of GSE251778"),
           size = 2.75, lineheight = 1.2, colour = ink) +
  annotate("text", x = 9.35, y = 3.55,
           label = "Promoter disease tests\nFull-CpG EWAS",
           fontface = "bold", size = 2.75, lineheight = 1.2, colour = slate) +
  annotate("segment", x = 3.59, xend = 4.00, y = 6.88, yend = 6.88,
           arrow = grid::arrow(length = grid::unit(0.11, "in"), type = "closed"),
           linewidth = 0.65, colour = slate) +
  annotate("segment", x = 7.29, xend = 7.70, y = 6.88, yend = 6.88,
           arrow = grid::arrow(length = grid::unit(0.11, "in"), type = "closed"),
           linewidth = 0.65, colour = slate) +
  annotate("rect", xmin = 0.35, xmax = 10.95, ymin = 0.55, ymax = 1.85,
           fill = "white", colour = "#D6DEE2", linewidth = 0.5) +
  annotate("text", x = 5.65, y = 1.54, label = "ADDITIONAL TISSUE AND GENETIC CONTEXT",
           fontface = "bold", size = 2.8, colour = ink) +
  annotate("text", x = 0.70, y = 1.05,
           label = "Brain tissue: 7 regional comparisons across GSE53987, GSE54563 and GSE80655",
           hjust = 0, size = 2.25, colour = slate) +
  annotate("text", x = 0.70, y = 0.67,
           label = "Genetic context: MDD-anxiety LDSC and 25 targeted OXPHOS MAGMA sets",
           hjust = 0, size = 2.25, colour = slate) +
  annotate("text", x = 5.65, y = 0.15,
           label = "Arrows show the exploratory analysis sequence, not a causal pathway.",
           size = 2.55, colour = slate, fontface = "italic") +
  coord_cartesian(xlim = c(0, 11.3), ylim = c(-0.1, 7.2), expand = FALSE, clip = "off") +
  theme_void(base_family = "sans") +
  theme(plot.margin = margin(8, 10, 14, 10))
save_figure(fig1, "Figure1_study_progression", 7.2, 4.7)

# Figure 2: all eight cohort-native OXPHOS rows from the existing GO-BP screen.
nes <- read_tsv(file.path(native, "Figure_2_native_GO_BP_NES.tsv"))
stopifnot(nrow(nes) == 8, all(nes$ID == "GO:0006119"))
nes$Display_label[is.na(nes$Display_label) | nes$Display_label == ""] <- nes$Analysis[is.na(nes$Display_label) | nes$Display_label == ""]
nes$Display_label <- factor(nes$Display_label, levels = rev(unique(nes$Display_label)))
nes$q_text <- paste0("q = ", vapply(nes$p.adjust, q_label, character(1)),
                     "   |   genes = ", nes$setSize)
fig2 <- ggplot(nes, aes(x = NES, y = Display_label)) +
  geom_vline(xintercept = 0, linewidth = 0.45, colour = "#8B979E") +
  geom_segment(aes(x = 0, xend = NES, yend = Display_label),
               linewidth = 0.85, colour = "#D4DDE1") +
  geom_point(aes(colour = ifelse(NES < 0, "Down", "Up"),
                 shape = ifelse(grepl("female|male", Display_label, ignore.case = TRUE),
                               "Stratified", "Overall")),
             size = 2.8, stroke = 0.85) +
  geom_text(aes(x = 3.1, label = q_text), hjust = 0, size = 2.65, colour = ink) +
  scale_colour_manual(values = c(Down = blue, Up = copper)) +
  scale_shape_manual(values = c(Overall = 16, Stratified = 15)) +
  scale_x_continuous(breaks = seq(-2, 2, 1), limits = c(-3.0, 5.1),
                     expand = expansion(mult = c(0.01, 0.01))) +
  labs(title = "Exploratory screening nominated OXPHOS",
       subtitle = "Cohort-native GO-BP NES; positive values indicate upward MDD enrichment",
       x = "Normalized enrichment score (NES)", y = NULL,
       colour = NULL, shape = NULL) +
  theme_paper(10) +
  theme(axis.text.y = element_text(size = 8.7),
        legend.position = "top",
        plot.margin = margin(8, 16, 8, 8)) +
  guides(colour = guide_legend(order = 1), shape = guide_legend(order = 2)) +
  coord_cartesian(clip = "off")
save_figure(fig2, "Figure2_exploratory_OXPHOS_NES", 7.2, 4.6)

# Figure 3: exact paired CAMERA P values, with FRY in a separate panel.
correlation_order <- read_tsv(file.path(tab, "Figure2_display_order.tsv"))
corr <- read_tsv(file.path(tab, "Table2_correlation_comparisons_11.tsv"))
stopifnot(nrow(corr) == 11)
corr$Display_order <- correlation_order$Display_order[match(corr$Comparison,
                                                              correlation_order$Comparison)]
stopifnot(!anyNA(corr$Display_order))
corr <- corr[order(corr$Display_order), ]
corr$rho_label <- paste0(corr$Comparison, "  (rho = ", sprintf("%.2f", corr$Correlation), ")")
corr$rho_label <- factor(corr$rho_label, levels = rev(corr$rho_label))
corr$y <- as.numeric(corr$rho_label)
corr$Fixed_logP <- -log10(corr$Fixed_P)
corr$Estimated_logP <- -log10(corr$Estimated_P)
corr$FRY_logP <- -log10(corr$FRY_P)
corr_long <- rbind(
  data.frame(corr, Method = "Fixed correlation 0.01", logP = corr$Fixed_logP),
  data.frame(corr, Method = "Estimated correlation", logP = corr$Estimated_logP)
)
corr_long$Method <- factor(corr_long$Method, levels = c("Fixed correlation 0.01", "Estimated correlation"))
nominal <- -log10(0.05)
p3a <- ggplot(corr_long, aes(x = logP, y = rho_label)) +
  geom_vline(xintercept = nominal, linetype = "22", linewidth = 0.45, colour = "#8C969C") +
  geom_segment(data = corr, aes(x = Fixed_logP, xend = Estimated_logP,
                                y = rho_label, yend = rho_label),
               inherit.aes = FALSE, linewidth = 0.55, colour = "#CFD7DB") +
  geom_point(aes(shape = Method, fill = Method), size = 2.7, stroke = 0.8,
             colour = c(copper, blue)[as.integer(corr_long$Method)]) +
  scale_shape_manual(values = c("Fixed correlation 0.01" = 1, "Estimated correlation" = 21)) +
  scale_fill_manual(values = c("Fixed correlation 0.01" = "white", "Estimated correlation" = blue)) +
  scale_x_continuous(limits = c(0, 15), breaks = c(0, 1.3, 5, 10, 14),
                     labels = c("0", "1.3", "5", "10", "14"),
                     expand = expansion(mult = c(0, 0.02))) +
  labs(title = "CAMERA P values depend on the correlation setting",
       subtitle = "Paired tests retain the same cohort, expression input, set, model and contrast",
       x = expression(-log[10](P)), y = NULL) +
  theme_paper(9.5) +
  theme(axis.text.y = element_text(size = 7.7),
        legend.position = "top",
        plot.margin = margin(4, 6, 2, 8))
p3b <- ggplot(corr, aes(x = FRY_logP, y = rho_label)) +
  geom_vline(xintercept = nominal, linetype = "22", linewidth = 0.45, colour = "#8C969C") +
  geom_point(shape = 23, size = 2.6, fill = slate, colour = slate) +
  scale_x_continuous(limits = c(0, 1.75), breaks = c(0, nominal, 1.5),
                     labels = c("0", "1.3", "1.5"),
                     expand = expansion(mult = c(0, 0.02))) +
  labs(title = "Directional FRY", x = expression(-log[10](P)), y = NULL) +
  theme_paper(9.2) +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        plot.title = element_text(size = 9.3, face = "bold", colour = ink),
        plot.margin = margin(0, 6, 8, 8))
fig3 <- p3a / p3b + plot_layout(heights = c(3.0, 1.35), guides = "collect") +
  plot_annotation(tag_levels = "A") & theme(legend.position = "top")
save_figure(fig3, "Figure3_correlation_sensitivity", 7.2, 8.0)

# Figure 4: distinct participant-level score estimands remain in separate panels.
score <- read_tsv(file.path(tab, "Figure3_historical_score_display.tsv"))
score_order <- score$Label
score$Label <- factor(score$Label, levels = rev(score_order))
score$N_label <- c("n=192", "n=45", "n=67", "n=169", "k=3", "k=4")
score$Kind <- ifelse(grepl("synthesis|sensitivity", score$Label, ignore.case = TRUE),
                     "REML synthesis", "Cohort estimate")
score$y <- as.numeric(score$Label)
p4a <- ggplot(score, aes(y = y, x = Beta)) +
  geom_vline(xintercept = 0, colour = "#87949B", linewidth = 0.45) +
  geom_segment(aes(x = Low, xend = High, yend = y), linewidth = 0.72,
               colour = ifelse(score$Kind == "REML synthesis", ink, blue)) +
  geom_point(aes(shape = Kind, colour = Kind), size = 2.8, stroke = 0.8) +
  geom_text(aes(x = 0.53, label = N_label), hjust = 0, size = 2.7, colour = slate) +
  scale_shape_manual(values = c("Cohort estimate" = 16, "REML synthesis" = 18)) +
  scale_colour_manual(values = c("Cohort estimate" = blue, "REML synthesis" = ink)) +
  scale_y_continuous(breaks = seq_along(score_order),
                     labels = rev(score_order), expand = expansion(add = 0.55)) +
  scale_x_continuous(limits = c(-1.2, 0.86), breaks = seq(-1.0, 0.5, 0.5)) +
  labs(title = "A  Participant-level rank scores",
       subtitle = "Cohort estimates use HC3 normal intervals; syntheses use modified Hartung-Knapp intervals",
       x = "MDD-control coefficient (within-cohort SD units)", y = NULL) +
  theme_paper(9.2) +
  theme(axis.text.y = element_text(size = 8.2), legend.position = "top",
        plot.margin = margin(4, 28, 1, 8)) +
  coord_cartesian(clip = "off")

cell <- read_tsv(file.path(tab, "OXPHOS_cell_proxy_sensitivity.tsv"))
cell <- cell[as.character(cell$Primary_run) == "TRUE" &
               cell$Stratum == "overall" &
               cell$Signature == "OXPHOS_GO0006119" &
               cell$Model %in% c("base_covariates", "marker_proxy_adjusted"), ]
cohort_order <- c("GSE98793", "GSE39653", "GSE19738", "GSE251778")
cell$Cohort <- factor(cell$Cohort, levels = cohort_order)
cell$yy <- 5 - as.numeric(cell$Cohort)
cell$offset <- ifelse(cell$Model == "base_covariates", 0.12, -0.12)
cell$yy <- cell$yy + cell$offset
cell$Model_label <- ifelse(cell$Model == "base_covariates", "Base covariates", "Marker-proxy adjusted")
base_cell <- cell[cell$Model == "base_covariates", ]
adj_cell <- cell[cell$Model == "marker_proxy_adjusted", ]
joined_cell <- merge(base_cell[, c("Cohort", "Beta")], adj_cell[, c("Cohort", "Beta")],
                     by = "Cohort", suffixes = c("_base", "_adj"))
joined_cell$y <- 5 - as.numeric(factor(joined_cell$Cohort, levels = cohort_order))
p4b <- ggplot(cell, aes(x = Beta, y = yy)) +
  geom_vline(xintercept = 0, colour = "#87949B", linewidth = 0.45) +
  geom_segment(data = joined_cell, aes(x = Beta_base, xend = Beta_adj,
                                        y = y + 0.12, yend = y - 0.12),
               inherit.aes = FALSE, linewidth = 0.55, colour = "#D1D9DD") +
  geom_segment(aes(x = CI_low, xend = CI_high, yend = yy), linewidth = 0.65,
               colour = ifelse(cell$Model == "base_covariates", copper, blue)) +
  geom_point(aes(shape = Model_label, fill = Model_label), size = 2.6,
             colour = ink, stroke = 0.7) +
  scale_shape_manual(values = c("Base covariates" = 21, "Marker-proxy adjusted" = 21)) +
  scale_fill_manual(values = c("Base covariates" = "white", "Marker-proxy adjusted" = blue)) +
  scale_y_continuous(breaks = 4:1, labels = cohort_order, limits = c(0.5, 4.5)) +
  scale_x_continuous(limits = c(-0.48, 0.38), breaks = seq(-0.4, 0.2, 0.2)) +
  labs(title = "B  Expression-derived cell-marker sensitivity",
       subtitle = "Mean-z score coefficients; 95% HC3 t intervals",
       x = "MDD-control coefficient (standardized score units)", y = NULL) +
  theme_paper(9.2) +
  theme(axis.text.y = element_text(size = 8.2), legend.position = "top",
        plot.margin = margin(2, 8, 1, 8))

paired <- read_tsv(file.path(tab, "paired_RNA.tsv"))
paired <- paired[paired$Score == "RNA_OXPHOS_all", ]
stopifnot(nrow(paired) == 2)
paired$Model_label <- ifelse(paired$Model == "A_no_cell", "No cell adjustment", "EpiDISH-adjusted")
paired$y <- ifelse(paired$Model == "A_no_cell", 2, 1)
p4c <- ggplot(paired, aes(x = Beta, y = y)) +
  geom_vline(xintercept = 0, colour = "#87949B", linewidth = 0.45) +
  geom_segment(aes(x = CI_HC3_low, xend = CI_HC3_high, yend = y),
               linewidth = 0.65, colour = blue) +
  geom_point(shape = 21, size = 2.7, fill = blue, colour = ink, stroke = 0.65) +
  scale_y_continuous(breaks = c(2, 1),
                     labels = c("No cell adjustment", "EpiDISH-adjusted"),
                     limits = c(0.5, 2.5)) +
  scale_x_continuous(limits = c(-0.45, 0.55), breaks = seq(-0.4, 0.4, 0.2)) +
  labs(title = "C  Paired-subset RNA score",
       subtitle = "GSE251780 participants; 95% HC3 t intervals",
       x = "MDD-control coefficient (standardized score units)", y = NULL) +
  theme_paper(9.2) +
  theme(axis.text.y = element_text(size = 8.1), legend.position = "none",
        plot.margin = margin(2, 8, 8, 8))
fig4 <- p4a / p4b / p4c + plot_layout(heights = c(1.55, 1.35, 0.95))
save_figure(fig4, "Figure4_participant_score_estimates", 7.2, 8.3)

# Figure 5: only gene-region-paired, gene-averaged promoter disease coefficients.
prom <- read_tsv(file.path(tab, "promoter_DNAm.tsv"))
coverage <- read_tsv(file.path(tab, "promoter_coverage.tsv"))
prom <- merge(prom, coverage, by = "Set", all.x = TRUE, sort = FALSE)
set_order <- c("GO0006119_original_ID_repaired", "MitoCarta_complex_subunits_nuclear",
               "Complex_I", "Complex_II", "Complex_III", "Complex_IV", "Complex_V")
prom <- prom[prom$Set %in% set_order, ]
stopifnot(nrow(prom) == 14, all(!is.na(prom$Defined_genes)))
prom$Set <- factor(prom$Set, levels = set_order)
prom$set_y <- 8 - as.numeric(prom$Set)
prom$Model_label <- ifelse(prom$Model == "sex_adjusted", "Sex-adjusted", "Cell-adjusted")
prom$offset <- ifelse(prom$Model == "sex_adjusted", 0.12, -0.12)
prom$y <- prom$set_y + prom$offset
prom_sex <- prom[prom$Model == "sex_adjusted", c("Set", "Beta_score_SD", "set_y")]
prom_cell <- prom[prom$Model == "cell_adjusted", c("Set", "Beta_score_SD", "set_y")]
prom_pair <- merge(prom_sex, prom_cell, by = c("Set", "set_y"),
                   suffixes = c("_sex", "_cell"), sort = FALSE)
fig5 <- ggplot(prom, aes(x = Beta_score_SD, y = y)) +
  geom_vline(xintercept = 0, colour = "#87949B", linewidth = 0.45) +
  geom_segment(data = prom_pair,
               aes(x = Beta_score_SD_sex, xend = Beta_score_SD_cell,
                   y = set_y + 0.12, yend = set_y - 0.12),
               inherit.aes = FALSE, linewidth = 0.52, colour = "#D1D9DD") +
  geom_segment(aes(x = CI_low, xend = CI_high, yend = y), linewidth = 0.67,
               colour = ifelse(prom$Model == "sex_adjusted", copper, blue)) +
  geom_point(aes(shape = Model_label, fill = Model_label), size = 2.65,
             colour = ink, stroke = 0.72) +
  scale_shape_manual(values = c("Sex-adjusted" = 21, "Cell-adjusted" = 21)) +
  scale_fill_manual(values = c("Sex-adjusted" = "white", "Cell-adjusted" = blue)) +
  scale_y_continuous(breaks = 7:1,
                     labels = c("GO:0006119 OXPHOS", "Nuclear MitoCarta subunits",
                                "Complex I", "Complex II", "Complex III", "Complex IV", "Complex V"),
                     limits = c(0.5, 7.5)) +
  scale_x_continuous(limits = c(-0.43, 0.48), breaks = seq(-0.4, 0.4, 0.2)) +
  labs(title = "Promoter disease associations were imprecise across definitions",
       subtitle = "Gene-region-paired, gene-averaged scores in 147 participants; all intervals are 95% HC3",
       x = "MDD coefficient (standardized promoter-score units)", y = NULL) +
  theme_paper(9.6) +
  theme(axis.text.y = element_text(size = 8.8), legend.position = "top",
        plot.margin = margin(8, 12, 8, 8))
save_figure(fig5, "Figure5_paired_promoter_disease_associations", 7.2, 5.6)

# Copy all numerical source tables used in the plots without modifying their values.
copy_source(file.path(native, "Figure_2_native_GO_BP_NES.tsv"), "Figure2_native_GO_BP_NES.tsv")
copy_source(file.path(tab, "Table2_correlation_comparisons_11.tsv"), "Figure3_correlation_comparisons_11.tsv")
copy_source(file.path(tab, "Figure3_historical_score_display.tsv"), "Figure4_rank_score_display.tsv")
copy_source(file.path(tab, "OXPHOS_cell_proxy_sensitivity.tsv"), "Figure4_cell_proxy_sensitivity.tsv")
copy_source(file.path(tab, "paired_RNA.tsv"), "Figure4_paired_RNA_scores.tsv")
copy_source(file.path(tab, "promoter_DNAm.tsv"), "Figure5_promoter_disease_models.tsv")
copy_source(file.path(tab, "promoter_coverage.tsv"), "Figure5_promoter_coverage.tsv")

inventory <- data.frame(
  Figure = 1:5,
  Main_panel = c("Study progression and evidence layers",
                 "Exploratory cohort-native OXPHOS NES",
                 "Fixed/estimated CAMERA and directional FRY",
                 "Rank scores, cell-marker sensitivity, paired-subset RNA score",
                 "Gene-region-paired promoter disease scores"),
  Source = c("Table 1 extracted from v6.12 manuscript",
             "Figure_2_native_GO_BP_NES.tsv",
             "Table2_correlation_comparisons_11.tsv",
             "Figure3_historical_score_display.tsv; OXPHOS_cell_proxy_sensitivity.tsv; paired_RNA.tsv",
             "promoter_DNAm.tsv; promoter_coverage.tsv"),
  Transformation = c("Layout only; cohort counts and labels read from Table 1",
                     "Plot existing NES and within-run BH q values; no recalculation",
                     "Display -log10 of saved P values; FRY remains a separate panel",
                     "Display saved coefficients and confidence intervals in separate panels",
                     "Display saved coefficients and HC3 intervals for the repaired paired definition"),
  stringsAsFactors = FALSE
)
write.table(inventory, file.path(out_dir, "FIGURE_SOURCE_INVENTORY_v6_12.tsv"),
            sep = "\t", row.names = FALSE, quote = TRUE)

writeLines(capture.output(sessionInfo()), file.path(out_dir, "figure_render_sessionInfo.txt"))
cat("Rendered Figures 1–5 from existing v6.11 result tables. No model fitting was performed.\n")
