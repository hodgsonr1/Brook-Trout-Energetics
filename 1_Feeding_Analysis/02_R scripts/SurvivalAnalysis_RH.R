###########
## Survival Analysis and Learning ##

#Purpose: 
  #Code to conduct survival analysis on data from Kenauk 2024 and 2025 latency to feed experiments at 15C and 10C 
  #following various simulated C&R events on brook trout (salvenlinus frontinalis)

##Author: Ryan Hodgson
## Date created: Nov 22 2024
## Date last modified: Aug 21 2026
####################

#load packages
library(survival)
library(tidyverse)
library(lubridate)
library(survminer)
library(gtsummary)
library(flextable) 
library(car)
library(purrr)
library(coxme)

#load data
feed_df <- read.csv("01_Feeding_Exp_Data/FeedingExpTrialData_Kenauk_Combined_14_07_2025.csv", header=TRUE)

### CLEAN DATA -----------
#create a column for feed time in hours 
feed_df$Latency_hrs<- feed_df$Feeding.latency..mins./60

#change temp and treatments to factor 
feed_df$Temperature..C <- factor(feed_df$Temperature..C)
feed_df$Treatment   <- factor(feed_df$Treatment)

#make weight on a meaningful scale (kg)
feed_df$Weight.kg <- feed_df$Weight..g./1000

#drop rows with weight NA
feed_df <- feed_df[!is.na(feed_df$Weight.kg), ]

###################################################
#### Descriptive Statistics ######
###################################################
#overall weight and Tl summary
feed <- feed_df %>%
  summarise(
  n=n(),
  mean_weight = mean(Weight..g.),
  sd_weight = sd(Weight..g.),
  mean_tl = mean(Total.Length..mm.),
  sd_tl = sd(Total.Length..mm.))
print(feed)

#summary stats for each group
feed_summary <- feed_df %>%
  group_by(Temperature..C, Treatment) %>%
  summarise(
    n = n(),
    mean_latency = mean(Latency_hrs, na.rm = TRUE),
    sd_latency = sd(Latency_hrs, na.rm = TRUE),
    se_latency = sd_latency / sqrt(n)
  ) %>%
  arrange(Temperature..C, Treatment)
feed_summary
#summary stats for 15C fish
summary15 <- feed_df %>%
  filter(Temperature..C == 15) %>%
  summarise(
    n = n(),
    mean_latency = mean(Latency_hrs, na.rm = TRUE),
    sd_latency = sd(Latency_hrs, na.rm = TRUE),
    se_latency = sd_latency / sqrt(n)
  )
summary15
#summary stats for 10C fish
summary10 <- feed_df %>%
  filter(Temperature..C == 10) %>%
  summarise(
    n = n(),
    mean_latency = mean(Latency_hrs, na.rm = TRUE),
    sd_latency = sd(Latency_hrs, na.rm = TRUE),
    se_latency = sd_latency / sqrt(n)
  )
summary10

#body weight checks 
glimpse(feed_df)
## boxplot for body mass
ggplot(feed_df, aes(x = Treatment, y = Weight..g., fill = (Temperature..C))) +
  geom_boxplot(position = position_dodge(width = 0.8)) +
  labs(
    title = "Boxplot of Weight by Treatment and TempMeta",
    x = "Treatment",
    y = "Weight (g)",
    fill = "TempMeta"
  ) +
  theme_minimal()
### test to make sure no difference in size across groups 
leveneTest(Weight..g. ~ Treatment * Temperature..C, data = feed_df)# unequal variances!
#equal variance
#ANOVA
weight <- aov(Weight..g. ~ Temperature..C * Treatment, data = feed_df)
summary(weight)
tl <-aov(Total.Length..mm. ~ Temperature..C * Treatment, data = feed_df)
summary(tl)

#no significant differences in weight or Tl across treatments

###################################################
#### COX - PROPORTIONAL HAZARD SURV ANALYSIS ######
###################################################

#order predictors correctly
feed_df$Temperature..C <- relevel(feed_df$Temperature..C, ref = "10")
feed_df$Treatment <- relevel(feed_df$Treatment, ref = "control_A")

#Fit mixed effects cox model 
m_full <- coxme(
  Surv(Latency_hrs, Feeding_Event.Y.N.) ~ 
    (Treatment) + (Temperature..C) + (Weight.kg) + (1 | Trial..),#include random effect of trial
  data = feed_df 
)
summary(m_full)

##Omnibus Likelihood ratio test
#weight
m_no_weight <- coxme(
  Surv(Latency_hrs, Feeding_Event.Y.N.) ~ Treatment + Temperature..C + (1 | Trial..),
  data = feed_df
)
anova(m_no_weight, m_full)
#temperature
m_no_temp <- coxme(
  Surv(Latency_hrs, Feeding_Event.Y.N.) ~ Treatment + Weight.kg + (1 | Trial..),
  data = feed_df
)
anova(m_no_temp, m_full)
#treatment
m_no_treat <- coxme(
  Surv(Latency_hrs, Feeding_Event.Y.N.) ~ Temperature..C + Weight.kg + (1 | Trial..),
  data = feed_df
)
anova(m_no_treat, m_full)



#------------------------------------------------------
#check assumptions ######
# 1) Proportional hazards assumption
test.ph <- cox.zph(m_full)
test.ph
ggcoxzph(test.ph) 
#schoenfield residuals look relatively flat and are within bounds of +-2 SE from fit. 
#visually no systematic pattern of effect of weight over time.


# 2) Linearity (Weight.kg) 
fit_ph <- coxph(Surv(Latency_hrs, Feeding_Event.Y.N.) ~ Weight.kg + cluster(Trial..),
                data = feed_df)

Y <- residuals(fit_ph, type = "martingale")
X <- feed_df$Weight.kg

plot(X, Y, pch = 20, col = "darkgray",
     ylab = "Martingale residual",
     xlab = "Weight (kg)",
     main = "Martingale residuals vs Weight")
abline(h = 0, lty = 3)
lines(smooth.spline(X, Y, df = 7), lty = 2, lwd = 2)



##############################################################################
#Calculating feeding inputs for FB4 model
##############################################################################
# this code finds the fraction of time feeding over 24hrs after stress from C&R/Temp.
#description of steps:
#(1) Fit Cox model without random effects; use it to generate covariate-specific predicted survival curves.

#(2) Create function to calculate F_after for a given covariate profile
#    (temperature, treatment, weight) from its predicted survival curve
#    F_after represents the fraction of 24hrs spent feeding during recovery from C&R.

#(3) Calculate F_after for each temperature group (10 & 15C) based on a median fish (~300g) weight.

#(4) Quantify uncertainty of F_after values 
  #Conduct a non-parametric bootstrap of Cox Model: 
  #re sample data with replacement (same n as original) 1000 times
  #refit the cox model on each re-sample. 
  #from each refit, compute F_after at both temperatures (10 &15C) with fixed fish weight (300g)

#(5) Save results to FB4 model for use in bioenergetics calculations


#fit cox model without random effect 
#use to fit KM curves for weight x temp groupings and extract mean time to return to feeding for FB4

cph <- coxph(
  Surv(Latency_hrs, Feeding_Event.Y.N.) ~ 
    (Treatment) * (Temperature..C) + (Weight.kg),
  data = feed_df 
)
summary(cph)

# Function to calculate F_after (fraction of day feeding) for a given covariate profile.
calc_individual_p_stress <- function(temp, treat, weight, model) {
  nd <- tibble(
    Temperature..C = temp,
    Treatment = treat,
    Weight.kg = weight
  )

  # Fit survival curve for this individual
  sf <- survfit(model, newdata = nd)

  # Time grid for integration over 24 h
  times <- seq(0, 24, by = 0.25)

  # Extract S(t) for each time point
  S_t <- summary(sf, times = times)$surv

  # summary.survfit() truncates rather than extrapolating when the longest
  # observed time in the (possibly resampled) data falls short of 24h - pad
  # any missing trailing points flat at the last available value, otherwise
  # S_t silently comes back shorter than `times` and the integration below
  # would multiply mismatched-length vectors via R's recycling rule.
  if (length(S_t) < length(times)) {
    last_val <- if (length(S_t) > 0) S_t[length(S_t)] else 1
    S_t <- c(S_t, rep(last_val, length(times) - length(S_t)))
  }
  # Replace any NA (before first event) with 1
  S_t[is.na(S_t)] <- 1

  # Trapezoidal integration of S(t) over [0, 24hrs]
  dt <- diff(times)
  S_mid <- (S_t[-1] + S_t[-length(S_t)]) / 2
  int_S <- sum(S_mid * dt)

  # Fraction of the day feeding (1 - fraction not feeding)
  F_after <- 1 - (int_S / 24)

  return(F_after)
}

# Weight classes (Small/Median/Large) were dropped from the downstream growth model:
# Weight.kg is a large, significant predictor here (see summary(cph) above), but the
# effect is very likely competitive exclusion by larger tankmates in the 4-fish holding
# groups (not a physiological/behavioural stress recovery-rate difference a wild fish would show),
# and the downstream p0 calibration (Chalifoux, 113-299g) rarely reaches the old Large
# class (481.6g) either. F_after is now predicted at a single representative weight:
# the dataset's own median, which is also within Chalifoux's observed range.
median_weight_kg <- unname(quantile(feed_df$Weight.kg, probs = 0.50, na.rm = TRUE))
cat("\nRepresentative weight for F_after (dataset median, kg):", median_weight_kg, "\n")

# Calculate F_after for each temperature at the single representative weight
p_24h_tbl <- tibble(Temperature..C = c("10", "15")) %>%
  rowwise() %>%
  mutate(
    Weight_kg = median_weight_kg,
    F_after = calc_individual_p_stress(
      factor(Temperature..C, levels = levels(feed_df$Temperature..C)),
      factor("control_A", levels = levels(feed_df$Treatment)), 
      Weight_kg,
      cph
    )
  ) %>%
  ungroup() %>%
  select(Temperature..C, Weight_kg, F_after)

print(p_24h_tbl)

#save to disk for use in bioenergetics model
saveRDS(p_24h_tbl, "../3_BT_FB4_model/Inputs/p_24h_tbl.rds")

##############################################################################
# Bootstrap F_after: propagate Cox model parameter uncertainty ----
##############################################################################
# Resamples the raw trial rows (not the fitted model), refits the same coxph
# model each time, and recomputes F_after at the SAME fixed representative
# weight - only the feeding-resumption relationship's uncertainty propagates,

N_boot <- 1000
set.seed(42)

fit_boot_F_after <- function(resampled) {
  cph_b <- tryCatch(
    coxph(Surv(Latency_hrs, Feeding_Event.Y.N.) ~ Treatment * Temperature..C + Weight.kg,
          data = resampled),
    error = function(e) NULL, warning = function(w) NULL
  )
  if (is.null(cph_b) || any(is.na(coef(cph_b)))) return(NULL)

  tryCatch({
    tibble(Temperature..C = c("10", "15")) %>%
      rowwise() %>%
      mutate(
        F_after = calc_individual_p_stress(
          factor(Temperature..C, levels = levels(feed_df$Temperature..C)),
          factor("control_A", levels = levels(feed_df$Treatment)),
          median_weight_kg, cph_b
        )
      ) %>%
      ungroup()
  }, error = function(e) NULL)
}

boot_results <- map(seq_len(N_boot), function(i) {
  resampled <- feed_df[sample(nrow(feed_df), nrow(feed_df), replace = TRUE), ]
  # guard against a degenerate resample losing an entire factor level
  resampled$Temperature..C <- factor(resampled$Temperature..C, levels = levels(feed_df$Temperature..C))
  resampled$Treatment <- factor(resampled$Treatment, levels = levels(feed_df$Treatment))
  fit_boot_F_after(resampled)
})

n_failed <- sum(map_lgl(boot_results, is.null))
cat("\nBootstrap refits failed/skipped:", n_failed, "of", N_boot, "\n")

F_after_boot <- boot_results %>%
  compact() %>%
  imap_dfr(~ mutate(.x, boot_iter = .y))

cat("\n=== F_after bootstrap summary (n =", length(unique(F_after_boot$boot_iter)), "successful refits) ===\n")
print(F_after_boot %>% group_by(Temperature..C) %>%
        summarise(mean_F_after = mean(F_after), sd_F_after = sd(F_after),
                  lo = quantile(F_after, .025), hi = quantile(F_after, .975), .groups = "drop"))

saveRDS(F_after_boot, "../3_BT_FB4_model/Inputs/F_after_boot.rds")



