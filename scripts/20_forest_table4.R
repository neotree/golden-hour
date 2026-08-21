# =============================================================================
# Golden Hour Analysis -- SMCH Zimbabwe
# JULY 2026 REQUESTS -- Script 20: FOREST PLOT of the Table 4/5 odds ratios
# =============================================================================
# Publication-quality forest plot of the per-baby mortality odds ratios, read
# straight from Script 18's output (18a_or_full_models.csv). Shows all three
# strata (all eligible / LBW <2500 g / sensitivity 1500-2499 g), both exposure
# comparisons (DCC only, Both), and the three models (unadjusted, adjusted for
# birth weight, adjusted for birth weight + sex). Reference: neither DCC nor ESSC.
# ESSC-only is omitted (reported as counts only; see Table 3).
#
# INPUT : outputs/<DATA_TAG>/18a_or_full_models.csv  (run Script 18 first)
# OUTPUT: outputs/<DATA_TAG>/20_forest_table4.png
#
# DSH note: ASCII-only.
# =============================================================================

library(tidyverse)

args <- commandArgs(trailingOnly = FALSE)
sf <- args[grep("--file=", args)]
SCRIPT_DIR <- if (length(sf) > 0)
  dirname(normalizePath(sub("--file=", "", sf), mustWork = FALSE)) else getwd()

DATA_TAG <- "20260525_apgar1"   # which Script 18 run to plot (matches DATA_FILE date)
IN  <- file.path(SCRIPT_DIR, "outputs", DATA_TAG, "18a_or_full_models.csv")
OUT <- file.path(SCRIPT_DIR, "outputs", DATA_TAG, "20_forest_table4.png")

if (!file.exists(IN)) stop("Not found: ", IN, " -- run Script 18 first.")

d <- read_csv(IN, show_col_types = FALSE) %>%
  mutate(
    stratum = factor(stratum,
      levels = c("All eligible infants", "LBW <2500 g", "Sensitivity 1500-2499 g")),
    comparison = factor(comparison, levels = c("Both DCC and ESSC", "DCC only")),
    model = recode(model,
      "unadjusted"        = "Unadjusted",
      "adjusted (BW)"     = "Adjusted (BW)",
      "adjusted (BW+sex)" = "Adjusted (BW + sex)"),
    model = factor(model, levels = c("Unadjusted", "Adjusted (BW)", "Adjusted (BW + sex)")))

pal <- c("Unadjusted" = "#BBBBBB", "Adjusted (BW)" = "#4477AA",
         "Adjusted (BW + sex)" = "#228833")

p <- ggplot(d, aes(x = OR, y = comparison, colour = model)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey40") +
  geom_pointrange(aes(xmin = CI_lower, xmax = CI_upper),
                  position = position_dodge(width = 0.6), linewidth = 0.6, size = 0.45) +
  facet_grid(rows = vars(stratum), switch = "y") +
  scale_x_log10(breaks = c(0.25, 0.5, 1, 2), limits = c(0.13, 2.2)) +
  scale_colour_manual(values = pal) +
  guides(colour = guide_legend(reverse = TRUE)) +
  labs(x = "Odds ratio for death before discharge (95% CI, log scale)",
       y = NULL, colour = NULL,
       title = "Death before discharge by Golden Hour intervention receipt",
       subtitle = paste0("Reference: received neither DCC nor ESSC.\n",
                         "ESSC-only omitted from the models (reported as counts; Table 3).")) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom",
        strip.placement = "outside",
        strip.text.y.left = element_text(angle = 0, face = "bold"),
        strip.background = element_rect(fill = "grey95", colour = NA),
        panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold", size = 12))

ggsave(OUT, p, width = 9.5, height = 5.4, dpi = 300)
message("Saved: ", OUT)
