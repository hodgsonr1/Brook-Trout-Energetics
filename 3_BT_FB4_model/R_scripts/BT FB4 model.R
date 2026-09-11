#Brook Trout (BT) bioenergetics model
#Purpose: 
#----  this code is built off a LMB bioenergetics model from JWB. 
#----- Code is updated and changed to fit parameters to BT.
#Authors:Ryan Hodgson (RH) & Jake Brownscombe (JWB)
#last updated: Aug 21 2026
library(tidyverse)
library(patchwork)
library(flextable)
library(officer)

#basic model ----
oxycal <- 13560 #standard in FB4

# p0 (ration proportion), fish energy density (ED), and prey/feed energy density (EDP)
# are now calibrated from real Chalifoux 2025 hatchery data (see Chalifoux2025_Ref_Growth.R)
# rather than literature-range placeholders - both the Kenauk trial and the Chalifoux trial
# were captive brook trout fed commercial extruded salmonid pellets.
p0_boot <- readRDS("Inputs/p0_chalifoux_boot.rds")$p0   # 86 individual-fish p0 estimates
ED_boot <- readRDS("Inputs/ED_chalifoux_boot.rds")       # 6 tank final-day fish ED values (J/g)
EDP <- readRDS("Inputs/EDP_chalifoux.rds")               # representative pellet ED (J/g)

# F_after (post-capture feeding-suppression multiplier) bootstrap: 1000 refits of the
# Kenauk feeding-trial Cox model (see SurvivalAnalysis_RH.R), predicted at a single
# representative weight (0.3195 kg, the trial's own median) -
F_after_boot <- readRDS("Inputs/F_after_boot.rds")       # 1000 boot_iter x 2 temps
#####
param <- read.csv("Inputs/Parameters_official.csv")
bt <- param %>% filter(Species=="Brook Trout (juvenile & adult)")
bt$ED <- mean(ED_boot) #fallback scalar for the standalone diagnostic block below; MC loop draws its own ED per iteration
bt$ED <- as.numeric(bt$ED)

##############################------------
# Check and assess Brook Trout FB4 model components ------------
##############################------------

#data 
fb <- merge(data.frame(temp=seq(1,35,1)), data.frame(weight=seq(100,1000,5))) #create grid of possible temperatures and weights 
head(fb)

#metabolism 
#EQ2
Z <- log(bt$RQ) * (bt$RTM - bt$RTO)
Y <- log(bt$RQ) * (bt$RTM - bt$RTO + 2)
X <- (Z^2 * (1 + sqrt(1 + 40/Y))^2) / 400

V <- (bt$RTM - fb$temp) / (bt$RTM - bt$RTO)
fb$ft <- ifelse(fb$temp < bt$RTM, V^X * exp(X * (1 - V)), 1e-6)

Rmax <- bt$RA * fb$weight ^ bt$RB
fb$Met.J <- Rmax * fb$ft * fb$weight * oxycal

#metabolism values
ggplot(fb, aes(temp, ft))+geom_point()+
ggplot(fb, aes(temp, Rmax, col=weight))+geom_point()+
ggplot(fb, aes(temp, Met.J, col=weight))+geom_point()

#consumption
#EQ3 
CG1 <- (1/(bt$CTO-bt$CQ))*log((0.98*(1-bt$CK1))/(bt$CK1*0.02))
L1 <- exp(CG1*(fb$temp-bt$CQ))
KA <- (bt$CK1*L1) / (1 + bt$CK1*(L1-1))
CG2 <- (1/(bt$CTL-bt$CTM))*log((0.98*(1-bt$CK4))/(bt$CK4*0.02))
L2 <- exp(CG2*(bt$CTL-fb$temp))
KB <- (bt$CK4*L2) / (1 + bt$CK4*(L2-1))
fb$ft <- KA * KB #temp dependence curve
fb$Cmax <- bt$CA * (fb$weight) ^bt$CB #max consumption
fb$Cons.p <- 1 #max consumption at 1
fb$C <- fb$Cmax * fb$Cons.p * fb$ft #actual daily consumption.
fb$Cons.g <- fb$C*fb$weight #g of prey eaten
fb$Cons.J <- fb$Cons.g*EDP #joules

#consumption values:
ggplot(fb, aes(temp, ft))+geom_point()+
  ggplot(fb, aes(temp, Cmax, col=weight))+geom_point()+
  ggplot(fb, aes(temp, C, col=weight))+geom_point()+
  ggplot(fb, aes(temp, Cons.J, col=weight))+geom_point()

#wastes
#EG equation 2
# Egestion (fecal losses)
fb$Eg <- bt$FA * fb$temp^(bt$FB) * exp(bt$FG * fb$Cons.p) * fb$Cons.J
# Excretion (urinary losses)
fb$Ex <- bt$UA * fb$temp^(bt$UB) * exp(bt$UG * fb$Cons.p) * (fb$Cons.J - fb$Eg)

#SDA
fb$SDA <- bt$SDA *(fb$Cons.J-fb$Eg) 

#growth
fb$growth.J <- fb$Cons.J-fb$Met.J-fb$SDA-fb$Eg-fb$Ex
fb$growth.g <- fb$growth.J/bt$ED

#outputs
head(fb)
ggplot(fb, aes(temp, Cons.J, col=weight))+geom_point()+scale_color_viridis_c()+theme_bw()+
ggplot(fb, aes(temp, Met.J, col=weight))+geom_point()+scale_color_viridis_c()+theme_bw()+
ggplot(fb, aes(temp, SDA, col=weight))+geom_point()+scale_color_viridis_c()+theme_bw()+
ggplot(fb, aes(temp, Eg, col=weight))+geom_point()+scale_color_viridis_c()+theme_bw()+
ggplot(fb, aes(temp, Ex, col=weight))+geom_point()+scale_color_viridis_c()+theme_bw()+
ggplot(fb, aes(temp, growth.g, col=weight))+geom_point()+scale_color_viridis_c()+theme_bw()

##############################------------
# Running grow function ------------
##############################------------
#grow a fish over period of time and add simulated C&R 

#load data from resp and feeding projects
resp_df <- readRDS("Inputs/stress_EPOC.J.rds")
# F_after is dropped here - it's now drawn per Monte Carlo iteration from F_after_boot
# below (see the loop), not baked in as a single static value.
cons_df <- readRDS("Inputs/p_24h_tbl.rds") %>% select(-F_after)

##############################------------
# Set up Temp Df #
##############################------------
# Daily temperature regime: real Kenauk-property wild telemetry (Little Bent Lake,
# 8 control fish, see TemperatureInput.R) over 30-day window sits
# entirely within its well-sampled period (>=4 fish contributing on every day).
sim_template <- read.csv("Inputs/temp_daily_littlebent_30day.csv", header = TRUE)

##############################------------
#combine into 6 stress scenarios (3 EPOC x 2 temperatures; weight is now a single
# fixed representative value, not a scenario axis - see F_after_boot note above)
# Pivot EPOC data to long format
resp_long <- resp_df %>%
  pivot_longer(cols = everything(),
               names_to = "EPOC_level",
               values_to = "EPOC.J.g") %>%
  mutate(EPOC_level = case_when(
    EPOC_level == "EPOC.J.g_p10" ~ "low",
    EPOC_level == "EPOC.J.g_p50" ~ "med",
    EPOC_level == "EPOC.J.g_p90" ~ "high"
  ))

# Cross join to create all 6 combinations
stress_meta <- cons_df %>%
  crossing(resp_long)

# Fixed grid skeleton for all 24 scenarios (6 stress scenarios x 4 C&R event counts).
# No p0/p_stress/F_after columns here - p0, ED, and F_after are all Monte Carlo draws
# per iteration (p0 shared across the single representative weight, since it's a
# proportion of Cmax which already scales allometrically with weight), attached
# inside the MC loop below.
param_grid_skeleton <- stress_meta %>%
  crossing(n_events = c(0, 1, 3, 5))
head(param_grid_skeleton)

# Set up simulation parameters
ndays <- 30      # matches the Little Bent 30-day telemetry window (well-sampled, >=4 fish/day);
                 # still short of the 60-day Chalifoux calibration window, so p0 is not extrapolated
y_spacing <- 5   # days between C&R events - max n_events=5 schedules its last event on day 25,
                 # comfortably inside the 30-day window

# Load functions and run all 24 scenarios
source("R_scripts/functions.R")

##############################------------
# Monte Carlo propagation of p0, ED, and F_after uncertainty through all 24 scenarios ----
##############################------------
set.seed(42)
N_MC <- 1000  # final run

results <- map_dfr(seq_len(N_MC), function(iter) {
  p0_iter <- sample(p0_boot, 1)
  ED_iter <- sample(ED_boot, 1)
  meta_iter <- bt
  meta_iter$ED <- ED_iter

  # Draw one shared F_after bootstrap replicate and use both of its temperature
  # values together (not independent draws) - they come from a single Cox refit
  # and are correlated by construction.
  boot_idx <- sample(unique(F_after_boot$boot_iter), 1)
  F_after_iter <- F_after_boot %>%
    filter(boot_iter == boot_idx) %>%
    select(Temperature..C, F_after)

  grid_iter <- param_grid_skeleton %>%
    left_join(F_after_iter, by = "Temperature..C") %>%
    mutate(p_stress = F_after * p0_iter)  # informational only; run_scenario() recomputes this internally

  grid_iter %>%
    pmap_dfr(function(...) {
      scenario_row <- tibble(...)
      run_scenario(
        scenario_row = scenario_row,
        sim_template = sim_template,
        meta = meta_iter,
        ndays = ndays,
        p0 = p0_iter,
        y_spacing = y_spacing,
        oxycal = oxycal,
        EDP = EDP
      )
    }) %>%
    mutate(iter = iter, p0_used = p0_iter, ED_used = ED_iter, Fafter_boot_idx = boot_idx)
})

# View results
head(results)
glimpse(results)

#CALCULATE GROWTH DEVIATION + ENERGETIC COSTS (per iteration, so each iteration's
# angled scenarios are compared against that same iteration's own healthy/unstressed baseline)
outcomes <- results %>%
  group_by(iter) %>%
  mutate(
    # Baseline (healthy, n_events == 0) growth within this iteration
    baseline_growth_g = first(growth_g[n_events == 0]),
    baseline_energy    = first(Net_Energy[n_events == 0]),

    # Deviation from healthy baseline growth caused by the angling event(s), in grams
    growth_reduction_g = baseline_growth_g - growth_g,

    # Percent deviation from baseline growth (percent of healthy growth lost to angling)
    growth_reduction_pct = abs(growth_reduction_g / baseline_growth_g) * 100,

    # Energetic cost of the angling event(s): net energy lost vs. healthy baseline
    net_energy_cost = Net_Energy - baseline_energy,
  ) %>%
  group_by(iter, EPOC_level, Temperature) %>%
  mutate(
    # Net_Energy is only directly simulated by run_scenario() for n_events 0/1
    # (see functions.R), so net_energy_cost is NA for repeated captures (3, 5).
    # Approximate the total energetic cost as n_events x the single-capture cost:
    # each C&R event acts independently on its own scheduled day with no modeled
    # carryover/compounding between events, and weight/temperature are approximately
    # constant across the 30-day simulation window, so this is a reasonable
    # first-order scale-up consistent with how the model itself is structured.
    net_energy_cost = if_else(
      n_events %in% c(3, 5),
      net_energy_cost[n_events == 1] * n_events,
      net_energy_cost
    )
  ) %>%
  ungroup()

glimpse(outcomes)

##############################------------
# Aggregate across MC iterations: median + 95% CI per scenario combination ----
##############################------------
outcomes_summary <- outcomes %>%
  group_by(EPOC_level, Temperature, n_events) %>%
  summarise(
    # Raw percent growth (own-weight basis, not a baseline comparison) - retained for
    # the growth-by-n_events line plot further down
    pct_growth_med  = median(percent_growth), pct_growth_lo = quantile(percent_growth, .025), pct_growth_hi = quantile(percent_growth, .975),

    # Healthy (unstressed) fish weight gain, g
    healthy_growth_g_med = median(baseline_growth_g), healthy_growth_g_lo = quantile(baseline_growth_g, .025), healthy_growth_g_hi = quantile(baseline_growth_g, .975),
    # Angled (stressed) fish weight gain, g
    angled_growth_g_med  = median(growth_g), angled_growth_g_lo = quantile(growth_g, .025), angled_growth_g_hi = quantile(growth_g, .975),
    # Deviation from healthy baseline growth, g
    deviation_g_med      = median(growth_reduction_g), deviation_g_lo = quantile(growth_reduction_g, .025), deviation_g_hi = quantile(growth_reduction_g, .975),
    # Energetic cost of the angling event(s), J
    energy_cost_med      = median(net_energy_cost, na.rm = TRUE), energy_cost_lo = quantile(net_energy_cost, .025, na.rm = TRUE), energy_cost_hi = quantile(net_energy_cost, .975, na.rm = TRUE),
    # Percent deviation from baseline growth
    pct_deviation_med    = median(growth_reduction_pct, na.rm = TRUE), pct_deviation_lo = quantile(growth_reduction_pct, .025, na.rm = TRUE), pct_deviation_hi = quantile(growth_reduction_pct, .975, na.rm = TRUE),
    .groups = "drop"
  )

glimpse(outcomes_summary)

#SUMMARY TABLE DATA (median with 95% CI), across all C&R event counts
single_clean <- outcomes_summary %>%
  filter(n_events != 0) %>%
  mutate(
    `Deviation (g)`     = sprintf("%.1f (%.1f – %.1f)", deviation_g_med, deviation_g_lo, deviation_g_hi),
    `Percent Deviation` = sprintf("%.1f (%.1f – %.1f)", pct_deviation_med, pct_deviation_lo, pct_deviation_hi),
    `Energetic Cost`    = sprintf("%.1f (%.1f – %.1f)", energy_cost_med, energy_cost_lo, energy_cost_hi)
  ) %>%
  rename(`EPOC level` = EPOC_level, `Temperature C` = Temperature, `Number of Captures` = n_events) %>%
  arrange(`EPOC level`, `Temperature C`, `Number of Captures`) %>%
  select(
    `EPOC level`,
    `Temperature C`,
    `Number of Captures`,
    `Deviation (g)`,
    `Percent Deviation`,
    `Energetic Cost`
  )
ft <- flextable(single_clean)
ft <- autofit(ft)

doc <- read_docx()
doc <- body_add_flextable(doc, ft)
print(doc, target = "Graphs/revisions Sept 2026/single_captures_grouped.docx")

####### Growth deviation graph (median + 95% CI across MC iterations)
# EPOC_level is faceted explicitly here (not averaged away like the original point-estimate
# version) so the CI band reflects only MC sampling uncertainty (p0/ED/F_after), not EPOC
# blended in. Weight class is gone (single representative weight now) - Temperature moves
# to the x-axis/fill in its place.
plot <- ggplot(
  outcomes_summary %>% filter(n_events > 0),
  aes(x = Temperature,
      y = pct_deviation_med,
      fill = Temperature)
) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
  geom_errorbar(
    aes(ymin = pct_deviation_lo, ymax = pct_deviation_hi),
    position = position_dodge(width = 0.8), width = 0.3
  ) +
  facet_grid(
    rows = vars(EPOC_level),
    cols = vars(n_events),
    labeller = labeller(
      n_events   = function(x) paste(x, "C&R Events"),
      EPOC_level = function(x) paste("EPOC:", x)
    )
  ) +
  scale_fill_manual(
    values = c("10" = "steelblue", "15" = "coral"),
    name   = "Temperature (°C)"
  ) +
  labs(
    x     = "Recovery Temperature (°C)",
    y     = "Growth Deviation (%)"
  ) +
  theme_bw(base_size = 16) +
  theme(
    axis.text.x     = element_text(angle = 45, hjust = 1, face= "bold"),
    legend.position = "bottom",
    strip.background = element_rect(fill = "grey90", colour = "grey50"),
    strip.text      = element_text(face = "bold"),
    axis.text = element_text(face = "bold"),
    legend.text = element_text(face = "bold"),
    legend.title = element_text(face="bold"),
    axis.title = element_text(face="bold")
  )

plot
ggsave(
  filename = "Graphs/revisions Sept 2026/Pub/GrowthReduction_Percent.png",
  plot = plot,
  width = 7.5,      # inches (wider now that EPOC is a facet row too)
  height = 6,       # inches
  dpi = 300
)

### line plots for growth reduction w. baseline growth included.
#Growth by temperature x EPOC level x repeated captures (single representative weight)
# median line + 95% CI ribbon across MC iterations (p0/ED/F_after sampling uncertainty)
percent <- ggplot(
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

ggsave(
  filename = "Graphs/revisions Sept 2026/Pub/Percent_Growth.png",
  plot = percent,
  width = 6.5,      # inches
  height = 4.5,     # inches
  dpi = 300
)
print(percent)
