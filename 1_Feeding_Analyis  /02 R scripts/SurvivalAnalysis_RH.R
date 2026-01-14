###########
## Survival Analysis and Learning ##

#Purpose: 
  #Code to conduct survival analysis on data from Kenauk 2024 and 2025 latency to feed experiments at 15C and 10C 
  #following various simulated C&R events on brook trout (salvenlinus frontinalis)

#Instructions: 
  #1) open R project file "feeding exp"
  #2) read the following script, loading in the meta data file containing the latency to feed data obtained from scored experimental videos.
  #3) run code to conduct CPH analysis and check assumptions 

##Author: Ryan Hodgson
## Date created: Nov 22 2024
## Date last modified: Nov 9 2025
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
feed_df <- read.csv("01 Feeding_Exp_Data/FeedingExpTrialData_Kenauk_Combined_14_07_2025.csv", header=TRUE)

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
test.ph <- cox.zph(mixed)
test.ph
ggcoxzph(test.ph) 
#schoenfield residuals look flat and are within bounds of +-2 SE from fit. 
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
#(1) Calculate weight percentiles (10th, 50th, 90th) from the dataset
#(2) For each Temperature (10C, 15C), calculate F_after for each weight class (Small, Median, Large)
#(3) F_after represents the fraction of 24hrs spent feeding during recovery from C&R
#(4) Save results to FB4 model for use in bioenergetics calculations

#fit cox model without random effect 
#use to fit KM curves for weight x temp groupings and extract mean time to return to feeding for FB4

cph <- coxph(
  Surv(Latency_hrs, Feeding_Event.Y.N.) ~ 
    (Treatment) * (Temperature..C) + (Weight.kg),
  data = feed_df 
)
summary(cph)


### Calculate survival predictions for weight classes by temperature ###
# Use control_A as reference treatment (no significant treatment effect)
# Calculate F_after for small (10th percentile), median, and large (90th percentile) fish at each temperature

# Function to calculate F_after (fraction of day feeding) for a single individual
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

# Calculate weight percentiles from the dataset
weight_percentiles <- quantile(feed_df$Weight.kg, probs = c(0.10, 0.50, 0.90), na.rm = TRUE)

# Create table with weight classes
weight_classes <- tibble(
  Weight_Class = c("Small", "Median", "Large"),
  Weight_kg = c(weight_percentiles[1], weight_percentiles[2], weight_percentiles[3])
)

# Calculate F_after for each temperature and weight class combination
p_24h_tbl <- expand.grid(
  Temperature..C = c("10", "15"),
  Weight_Class = c("Small", "Median", "Large"),
  stringsAsFactors = FALSE
) %>%
  as_tibble() %>%
  left_join(weight_classes, by = "Weight_Class") %>%
  rowwise() %>%
  mutate(
    F_after = calc_individual_p_stress(
      factor(Temperature..C, levels = levels(feed_df$Temperature..C)),
      factor("control_A", levels = levels(feed_df$Treatment)),
      Weight_kg,
      cph
    )
  ) %>%
  ungroup() %>%
  select(Temperature..C, Weight_Class, Weight_kg, F_after)

print(p_24h_tbl)

#save to disk for use in bioenergetics model
saveRDS(p_24h_tbl, "../4_BT_FB4_model/Inputs/p_24h_tbl.rds")



