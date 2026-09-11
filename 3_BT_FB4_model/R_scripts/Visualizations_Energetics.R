## Visualizations - Brook Trout Energetics (FB4 Model) ##

# Description:
#   Visualizes outputs from Fish Bioenergetics 4 (FB4) simulations of catch-and-release (C&R)
#   scenarios for brook trout. Produces plots of growth reduction (%) across combinations
#   of recovery temperature, EPOC level, and number of C&R events, at a single representative
#   fish weight (weight classes were dropped - see BT FB4 model.R for rationale).
#   Uses `outcomes_summary` (median + 95% CI across Monte Carlo p0/ED/F_after draws), not the
#   raw per-iteration `outcomes` - since BT FB4 model.R now runs N_MC iterations, `outcomes`
#   has one row per iteration x scenario rather than one row per scenario.

# Instructions:
#   Run "BT FB4 model.R" first to generate the `outcomes_summary` data frame before running
#   this script.

# Author: Ryan Hodgson
# Date: Nov 29 2025
# Updated: Aug 20 2026 - dropped weight classes (single representative weight, 0.3195 kg) and
# switched to outcomes_summary (median + 95% CI) after BT FB4 model.R was changed to run
# p0/ED/F_after Monte Carlo iterations rather than a single point estimate.

# Dependencies:
#   tidyverse  - data wrangling (dplyr) and plotting (ggplot2)
#   flextable, officer - Net Energy Costs summary table

library(tidyverse)
library(flextable)
library(officer)

glimpse(outcomes_summary)
#Growth by temperature x EPOC level x repeated captures (single representative weight)
ggplot(
  outcomes_summary %>% arrange(n_events),
  aes(
    x = factor(n_events, levels = c(0, 1, 3, 5)),
    y = pct_growth_med,
    group = 1
  )
) +
  geom_ribbon(
    aes(ymin = pct_growth_lo, ymax = pct_growth_hi),
    alpha = 0.15, fill = "#5B8FD4", color = NA
  ) +
  geom_point(size = 3, color = "#5B8FD4") +
  geom_line(linewidth = 1, color = "#5B8FD4") +

  # Facets: temperature rows, EPOC columns
  facet_grid(
    rows = vars(Temperature),
    cols = vars(EPOC_level),
    labeller = labeller(
      Temperature = function(x) paste0(x, " °C"),
      EPOC_level  = function(x) paste("EPOC:", x)
    )
  ) +

  theme_bw(base_size = 14) +
  theme(
    strip.text        = element_text(size = 13, face = "bold"),
    strip.background  = element_rect(fill = "grey90", color = "black"),
    axis.title        = element_text(size = 14, face = "bold"),
    axis.text         = element_text(size = 12),
    panel.border      = element_rect(color = "black", linewidth = 0.8),
    panel.grid.minor  = element_blank(),
    panel.grid.major  = element_line(linewidth = 0.3, color = "grey85")
  ) +

  labs(
    x     = "Number of C&R Events",
    y     = "Percent Growth (%, median ± 95% CI)"
  )

#Faceted Barplots for Percent Deviation from Baseline Growth (median, with 95% CI error bars)
outcomes_filtered <- outcomes_summary %>% filter(n_events != 0)

ggplot(outcomes_filtered, aes(
  x = factor(n_events),
  y = pct_deviation_med,
  fill = Temperature
)) +
  geom_col(position = "dodge", alpha = 0.9) +
  geom_errorbar(
    aes(ymin = pct_deviation_lo, ymax = pct_deviation_hi),
    position = position_dodge(width = 0.9), width = 0.3
  ) +

  facet_grid(
    rows = vars(EPOC_level)
  ) +

  scale_fill_manual(
    values = c("10" = "steelblue", "15" = "coral"),
    name   = "Temperature (°C)"
  ) +

  theme_bw(base_size = 16) +
  theme(
    strip.text      = element_text(size = 18, face = "bold"),
    axis.title      = element_text(size = 20, face = "bold"),
    axis.text       = element_text(size = 16),
    axis.text.x     = element_text(angle = 45, hjust = 1),
    legend.title    = element_text(size = 18, face = "bold"),
    legend.text     = element_text(size = 16),

    panel.border    = element_rect(color = "black", linewidth = 1),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(linewidth = 0.4, linetype = "dashed")
  ) +

  labs(
    x = "C&R Events",
    y = "Per Cent Deviation from Baseline Growth (median, 95% CI)"
  )

#Summary Table: Net Energy Costs (all n_events, median with 95% CI)
# Net_Energy/net_energy_cost is only defined for n_events 0/1 by run_scenario()'s design
# (see functions.R), so 5/10-event rows are NA here by design, not a bug.
net_energy_table <- outcomes_summary %>%
  filter(!is.na(energy_cost_med)) %>%
  mutate(
    `Energy Cost J` = sprintf("%.1f (%.1f – %.1f)", energy_cost_med, energy_cost_lo, energy_cost_hi)
  ) %>%
  rename(EPOC = EPOC_level) %>%
  arrange(EPOC, Temperature, n_events) %>%
  select(Temperature, EPOC, n_events, `Energy Cost J`)

net_energy_ft <- flextable(net_energy_table) %>% autofit()

net_energy_doc <- read_docx()
net_energy_doc <- body_add_flextable(net_energy_doc, net_energy_ft)
print(net_energy_doc, target = "Graphs/revisions Sept 2026/NetEnergyCosts_summary.docx")

