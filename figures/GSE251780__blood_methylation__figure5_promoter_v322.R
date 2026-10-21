#!/usr/bin/env Rscript

# Re-export main-text Figure 5 from the formal v3.22.0 (v322) promoter results.
# Scope: Figure 5 only. Inputs are the frozen v322 tables (GO:0006119 = 148 defined genes).
# No models are refitted; this is display-only recomposition, mirroring the Figure 5
# main figure-export script's Figure 5 block but reading the v322 formal data files.
# PlanA (v3.21) promoter tables are retained unchanged as historical versions.

suppressPackageStartupMessages({
  library(ggplot2)
})

root <- "/Users/tan/Documents/Codex/2026-09-19/ple/outputs/manuscript_v6_12_optionB_2026-10-13"
out_dir <- file.path(root, "github_repository_draft")
fig_dir <- file.path(out_dir, "figures")
src_out <- file.path(out_dir, "figure_source_data")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(src_out, recursive = TRUE, showWarnings = FALSE)

# Formal v322 inputs (repository results directory)
tab <- file.path(root, "github_repository_draft/results/by_cohort/GSE251780")

ink <- "#25343D"; blue <- "#28789A"; copper <- "#B86D4B"
slate <- "#687984"

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

copy_source <- function(source, name) {
  target <- file.path(src_out, name)
  if (!file.copy(source, target, overwrite = TRUE)) stop("Could not copy source: ", source)
  target
}

# Figure 5: gene-region-paired, gene-averaged promoter disease coefficients (v322).
prom <- read_tsv(file.path(tab, "promoter_DNAm_v322.tsv"))
coverage <- read_tsv(file.path(tab, "promoter_coverage_v322.tsv"))
prom <- merge(prom, coverage, by = "Set", all.x = TRUE, sort = FALSE)
set_order <- c("GO0006119", "MitoCarta_complex_subunits_nuclear",
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

# Formal source data alongside the figure.
copy_source(file.path(tab, "promoter_DNAm_v322.tsv"), "Figure5_promoter_disease_models.tsv")
copy_source(file.path(tab, "promoter_coverage_v322.tsv"), "Figure5_promoter_coverage.tsv")

cat("Figure 5 (v322) exported to", fig_dir, "\n")
cat("Source data copied to", src_out, "\n")
print(prom[prom$Set == "GO0006119", c("Set", "Model", "Beta_score_SD", "CI_low", "CI_high", "P_HC3")])
