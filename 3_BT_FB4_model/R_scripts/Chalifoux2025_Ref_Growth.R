#Title 
### Chalifoux 2025 - Reference Growth 
#Script Author: Ryan Hodgson 

#Date Created: Aug 6 2026
#Date Last Modified: Aug 18 2026

#Description:
#using data from Chalifoux 2025, this script loads in baseline hatchery growth for control brook trout (salvelinus fontinalis)
# that are grown for 60 days under freshwater conditions.
# Replaces the VBGF-based reference growth (see vGBF_growth.R) as the target used to fit p0 in the BT FB4 model,
# since MNR long-term (BsM) field growth data is not a reliable analogue for captive/hatchery feeding conditions.

#Dependencies: (From Chalifoux 2025)
#Data files:
#1. R_wetmass.csv - Final and initial weights for BT over 60 days. Contains all data across all treatments from Chalifoux 2025
#2. Temperature_TX.csv` — daily water temperature  (for freshwater control tanks only)
#3. Ration_TX.csv` — daily grams of feed actually eaten (for freshwater control tanks only)
#4. Pred_E_TX.csv` — daily fish energy density trajectory (for freshwater control fish only)
#5. Prey_E.csv` — feed energy density (one file, all tanks shared the same feed)


#Data Citation:
#Chalifoux, Virginie, 2024, "Brook charr Wisconsin-type bioenergetics model evaluation", https://doi.org/10.5683/SP3/5A4LHI, Borealis, V2, UNF:6:XwzE3QVaC4di7r0G3W1sPg== [fileUNF]

library(tidyverse)
library(ggplot2)

getwd()

#Load Data ##-----------------------------
# semicolon-delimited; decimals use "." (not read.csv2 convention)
wet_mass <- read.csv("Inputs/Chalifoux_2025/R_wetmass.csv", header = TRUE, sep = ";")
glimpse(wet_mass)
view(wet_mass)
##############################------------
# Filter to freshwater (ED) control tanks ----
##############################------------
# traitement codes: ED = eau douce (freshwater control), ES = eau salee (saltwater treatment),
# INI = initial baseline sample (tank "0", measured before fish were allocated to tanks - not a real tank)
# -> ED tanks are the freshwater control group we want.

wet_mass_ED <- wet_mass %>%
  filter(traitement == "ED") %>%
  mutate(
    tank         = bassin,
    masse_init_g = masse_init_kg * 1000,
    masse_fin_g  = masse_fin_kg  * 1000
  ) %>%
  select(
    ID, tank,
    date_initiale, date_finale,
    masse_init_g, masse_fin_g,
    long_init_mm, long_fin_mm
  ) %>%
  arrange(tank, ID)

cat("\n=== Freshwater (ED) control tanks ===\n")
print(sort(unique(wet_mass_ED$tank)))

##############################------------
# Setup: bioenergetics parameters + model functions ----
##############################------------
oxycal <- 13560
param  <- read.csv("Inputs/Parameters_official.csv")
bt     <- param %>% filter(Species == "Brook Trout (juvenile & adult)")
source("R_scripts/functions.R")
# Note: bt$ED is left as the raw literature-range string here (unused) - every fit.p()
# call below passes an explicit per-day ED vector built from real Chalifoux measurements,
# so meta$ED is never actually read by grow() during calibration.

##############################------------
# Build per-tank daily profiles: temperature, prey energy density (EDP), fish ED ----
##############################------------
tanks <- c(3, 8, 10, 12, 14, 15)

prey_e <- read_csv("Inputs/Chalifoux_2025/Prey_E.csv", show_col_types = FALSE)
# columns: day, moulee (feed energy density, J/g) - step function, breakpoints at day 1/6/7/60

build_tank_profile <- function(tank_id) {
  temp_file <- sprintf("Inputs/Chalifoux_2025/Temperature/Temperature_T%d.csv", tank_id)
  pred_file <- sprintf("Inputs/Chalifoux_2025/Pred /Pred_E_T%d.csv", tank_id)  # NOTE: "Pred " folder has a trailing space

  temp_tbl <- read_csv(temp_file, show_col_types = FALSE) %>% arrange(day)
  ndays_tank <- max(temp_tbl$day)

  pred_tbl <- read_csv(pred_file, show_col_types = FALSE) %>% arrange(day)
  ed_start_day <- pred_tbl$day[1]; ed_start_val <- pred_tbl$predator[1]
  ed_end_day   <- pred_tbl$day[2]; ed_end_val   <- pred_tbl$predator[2]
  # ed_end_day is labeled 60 even for tanks (12/14/15) whose trial actually ran 61 days -
  # rule=2 extrapolates flat, so day 61 just gets the same ED as day 60.

  day_seq <- seq_len(ndays_tank)
  EDP_vec <- prey_e$moulee[findInterval(day_seq, prey_e$day)]
  ED_vec  <- approx(x = c(ed_start_day, ed_end_day), y = c(ed_start_val, ed_end_val),
                     xout = day_seq, rule = 2)$y

  tibble(
    tank = tank_id, day = day_seq,
    temp = temp_tbl$temperature[match(day_seq, temp_tbl$day)],
    EDP = EDP_vec, ED = ED_vec
  )
}

tank_profiles <- map_dfr(tanks, build_tank_profile)

tank_ED_final <- tank_profiles %>%
  group_by(tank) %>% slice_max(day, n = 1) %>% ungroup() %>%
  select(tank, ED_final = ED)
cat("\n=== Fish tissue ED, final day per tank (sanity check vs ~5390-5573) ===\n")
print(tank_ED_final)

##############################------------
# Fit p0 per individual fish (86 fish), using that fish's own tank's real profile ----
##############################------------
p0_fits <- wet_mass_ED %>%
  rowwise() %>%
  mutate(
    ndays_tank = max(tank_profiles$day[tank_profiles$tank == tank]),
    p0 = fit.p(
      initial_g    = masse_init_g,
      final_g      = masse_fin_g,
      sim_template = tank_profiles %>% filter(tank == .env$tank) %>% select(day, temp),
      meta         = bt,
      oxycal       = oxycal,
      EDP          = tank_profiles$EDP[tank_profiles$tank == tank],
      ndays        = ndays_tank,
      ED           = tank_profiles$ED[tank_profiles$tank == tank]
    )
  ) %>%
  ungroup()

##############################------------
# QC: convergence + p0 boundary checks (post-hoc, since fit.p() returns a bare scalar) ----
##############################------------
p0_fits <- p0_fits %>%
  rowwise() %>%
  mutate(
    check_final_w = {
      prof <- tank_profiles %>% filter(tank == .env$tank)
      d <- prof %>% select(day, temp)
      d$weight <- NA; d$E <- NA
      d$weight[1] <- masse_init_g
      d$E[1] <- masse_init_g * prof$ED[1]
      grow(
        data = d, meta = bt, ndays = ndays_tank, p0 = p0,
        stress_df = data.frame(day = integer(0), p_stress = numeric(0), EPOC.J.g = numeric(0)),
        oxycal = oxycal, EDP = prof$EDP, ED = prof$ED
      )[["final_weight"]]
    },
    converged = abs(check_final_w - masse_fin_g) <= 0.5,
    boundary_hit = p0 >= 0.995 | p0 <= 0.005
  ) %>%
  ungroup()

cat("\n=== p0 fit QC ===\n")
cat("Fish not converged (W.tol not met within max.iter):", sum(!p0_fits$converged), "of", nrow(p0_fits), "\n")
cat("Fish hitting p0 boundary [0,1] (model cannot explain observed growth within Cmax):",
    sum(p0_fits$boundary_hit), "\n")
print(p0_fits %>% filter(boundary_hit | !converged) %>%
        select(ID, tank, masse_init_g, masse_fin_g, p0, check_final_w))

##############################------------
# Ration-based validation: independent empirical p, cross-checked against bisection p0 ----
# (uses Ration_TX.csv; kept separate from the p0 Monte Carlo pool, not merged into it)
##############################------------
compute_ftC <- function(temp, meta) {
  CG1 <- (1 / (meta$CTO - meta$CQ)) * log((0.98 * (1 - meta$CK1)) / (meta$CK1 * 0.02))
  L1  <- exp(CG1 * (temp - meta$CQ))
  KA  <- (meta$CK1 * L1) / (1 + meta$CK1 * (L1 - 1))
  CG2 <- (1 / (meta$CTL - meta$CTM)) * log((0.98 * (1 - meta$CK4)) / (meta$CK4 * 0.02))
  L2  <- exp(CG2 * (meta$CTL - temp))
  KB  <- (meta$CK4 * L2) / (1 + meta$CK4 * (L2 - 1))
  KA * KB
}

ration_check <- map_dfr(tanks, function(tank_id) {
  ration <- read_csv(sprintf("Inputs/Chalifoux_2025/Ration/Ration_T%d.csv", tank_id), show_col_types = FALSE)
  prof   <- tank_profiles %>% filter(tank == tank_id)
  ndays_tank <- max(prof$day)

  tank_w <- wet_mass_ED %>% filter(tank == tank_id) %>%
    summarise(w0 = mean(masse_init_g), w1 = mean(masse_fin_g), n_fish = n())
  weight_traj <- approx(x = c(1, ndays_tank), y = c(tank_w$w0, tank_w$w1),
                         xout = seq_len(ndays_tank))$y

  d <- ration %>%
    filter(day <= ndays_tank) %>%
    left_join(prof %>% select(day, temp), by = "day") %>%
    mutate(
      weight = weight_traj[day],
      ftC = compute_ftC(temp, bt),
      Cmax_g = bt$CA * weight^bt$CB * ftC * weight,
      p_hat = ration / Cmax_g
    ) %>%
    filter(ration > 0, ration >= 0.9 * max(ration, na.rm = TRUE))  # data-driven plateau detector

  total_ration <- sum(ration$ration[ration$day <= ndays_tank], na.rm = TRUE)
  mass_gained  <- tank_w$w1 - tank_w$w0

  tibble(
    tank = tank_id, n_fish = tank_w$n_fish, n_plateau_days = nrow(d),
    p_hat_ration_mean = mean(d$p_hat), p_hat_ration_sd = sd(d$p_hat),
    total_ration_g = total_ration, mass_gained_g = mass_gained,
    FCR = total_ration / mass_gained
  )
})

cat("\n=== Per-tank FCR sanity check (total ration eaten / mass gained; expect ~1.0-2.5 for pellet feed) ===\n")
print(ration_check %>% select(tank, n_fish, total_ration_g, mass_gained_g, FCR))

bisection_by_tank <- p0_fits %>%
  group_by(tank) %>%
  summarise(p0_bisection_mean = mean(p0), p0_bisection_sd = sd(p0), .groups = "drop")

comparison <- left_join(ration_check, bisection_by_tank, by = "tank") %>%
  mutate(diff = p_hat_ration_mean - p0_bisection_mean)

cat("\n=== Bisection-based p0 vs Ration-based empirical p (validation cross-check) ===\n")
print(comparison %>% select(tank, p0_bisection_mean, p0_bisection_sd, p_hat_ration_mean, p_hat_ration_sd, diff))

##############################------------
# Save outputs for use by BT FB4 model.R ----
##############################------------
p0_chalifoux_boot <- p0_fits %>% ungroup() %>% select(ID, tank, p0)
saveRDS(p0_chalifoux_boot, "Inputs/p0_chalifoux_boot.rds")

ED_chalifoux_boot <- setNames(tank_ED_final$ED_final, tank_ED_final$tank)
saveRDS(ED_chalifoux_boot, "Inputs/ED_chalifoux_boot.rds")

# Representative EDP for the main Kenauk model: day-count-weighted mean of Prey_E's
# step-function values over a 60-day window (matching the Kenauk simulation length).
# The day-7 feed-switch timing is specific to Chalifoux's own trial protocol and has no
# meaningful mapping onto the Kenauk simulation calendar, so we use one fixed
# representative scalar rather than importing the step timing itself.
kenauk_ndays <- 60
EDP_daily_kenauk <- prey_e$moulee[findInterval(seq_len(kenauk_ndays), prey_e$day)]
EDP_chalifoux <- mean(EDP_daily_kenauk)
saveRDS(EDP_chalifoux, "Inputs/EDP_chalifoux.rds")
cat("\nRepresentative EDP for main model (Chalifoux pellet ED, 60-day weighted mean):", EDP_chalifoux, "J/g\n")

cat("\n=== p0 distribution summary (n =", nrow(p0_chalifoux_boot), ") ===\n")
print(summary(p0_chalifoux_boot$p0))
cat("SD:", sd(p0_chalifoux_boot$p0), "\n")

p0_hist <- ggplot(p0_chalifoux_boot, aes(p0)) +
  geom_histogram(bins = 20, fill = "steelblue", color = "white") +
  theme_bw(base_size = 14) +
  labs(x = "Fitted p0 (per fish, n=86)", y = "Count",
       title = "Chalifoux 2025 ED-tank individual p0 fits")
print(p0_hist)
ggsave("Graphs/revisions Sept 2026/p0_chalifoux_histogram.png", p0_hist, width = 6, height = 4.5, dpi = 300)









