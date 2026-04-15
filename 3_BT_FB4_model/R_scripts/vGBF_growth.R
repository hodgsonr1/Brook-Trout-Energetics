### Title: 
#MNR Growth Data Extraction ###

##Description: 
#load in brook trout age-growth data from MNR long term monitoring program. 
#use to build growth over time relationship using VGBF

#Ryan Hodgson 
#Nov 29 2025
#updated: Apr 6 2026

getwd()

#load packages
library(dplyr)
library(FSA)
library(FSAdata)
library(nlstools)
library(ggplot2)


#load data
growth <- read.csv("Inputs/BsM_BTdata.csv")

#clean data ----
# Select relevant columns and remove missing values
growth_clean <- growth %>%
  select(Age, RoundWeight, TotalLength) %>%
  filter(!is.na(Age) & !is.na(RoundWeight) & Age > 0 & RoundWeight > 0)
summary(growth_clean)

# Summary statistics
cat("Data Summary:\n")
cat("Number of observations:", nrow(growth_clean), "\n")
cat("Number of lakes:", length(unique(growth$Lake.Name)))
cat("Age range:", min(growth_clean$Age), "-", max(growth_clean$Age), "years\n")
cat("Weight range:", min(growth_clean$RoundWeight), "-", max(growth_clean$RoundWeight), "g\n")
cat("\nObservations per age class:\n")
print(table(growth_clean$Age))

#Estimate Starting Values 

#1) determine parameters from age-length data 
svTypical <- findGrowthStarts(TotalLength~Age,data=growth_clean)
print(svTypical)
Linf_S <- svTypical[["Linf"]]
K_S <- svTypical[["K"]]
t0_S <- svTypical[["t0"]]

#2) determine length-weight allometric relationship 
#regress the two, extract a and b. 
#coefficient (slope) = b
#intercept = a
lw_reg <- lm(log(RoundWeight)~log(TotalLength), data = growth_clean)
a <- exp(coef(lw_reg)[1])
b <- coef(lw_reg)[2]

#3) Determine Winf from a,b and Linf
Winf_S <- a * Linf_S^b
b_fixed <- b #fix allometric scaling coefficient to allow model convergence
print(Winf_S)
print(b_fixed)
# # Examine max values to check parameters
max_age <- max(growth_clean$Age)
max_weight <- max(growth_clean$RoundWeight)
mean_max_weight <- mean(growth_clean$RoundWeight[growth_clean$Age == max(growth_clean$Age)])
cat("\nMax observed age:", max_age, "years\n")
cat("Max observed weight:", max_weight, "g\n")
cat("Meam max observed weight:", mean_max_weight, "g\n")

#
cat("\nManual starting parameter values:\n")
cat("Winf:",Winf_S, "g\n")
cat("K:", K_S, "\n")
cat("t0:", t0_S, "\n")
cat("b", b, "\n")

#4) create list of starting paramters to fit model
start <- list(Winf = unname(Winf_S), K = unname(K_S), t0 = unname(t0_S))
#5) fit the VBGF using non-linear least squares
fit <- nls(RoundWeight ~ Winf * (1 - exp(-K * (Age - t0)))^b_fixed,
           data = growth_clean,
           start = start)

# Display model summary
cat("\n=== VBGF Model Summary ===\n")
summary(fit)

# Extract parameters
params <- coef(fit)
Winf <- params["Winf"]
K <- params["K"]
t0 <- params["t0"]
cat("\nFinal Model Parameter values:\n")
cat("Winf:",Winf, "g\n")
cat("K:", K, "\n")
cat("t0:", t0, "\n")
cat("b", b, "\n")

#check model fit 
par(mfrow = c(1, 2))
hist(residuals(fit), main = "Residuals", xlab = "Residuals")
plot(residuals(fit) ~ fitted(fit), main = "Residuals vs Fitted")
abline(h = 0, col = "red")
par(mfrow = c(1, 1))

#### ------ #
#generating confidence intervals using nlsBoot
bootTypical <- nlsBoot(fit, niter = 1000)
confint(bootTypical, plot = TRUE)
#
#VGBF visualization ----
# Create prediction data for smooth curve
age_pred <- seq(min(growth_clean$Age), max(growth_clean$Age), length.out = 100)
weight_pred <- predict(fit, newdata = data.frame(Age = age_pred))

pred_data <- data.frame(Age = age_pred, RoundWeight = weight_pred)

# Generate bootstrap CI band from nlsBoot parameter samples
boot_preds <- apply(bootTypical$coefboot, 1, function(p) {
  p["Winf"] * (1 - exp(-p["K"] * (age_pred - p["t0"])))^b_fixed
})
# boot_preds is a matrix: rows = age points, cols = bootstrap iterations
pred_data$lower <- apply(boot_preds, 1, quantile, probs = 0.025, na.rm = TRUE)
pred_data$upper <- apply(boot_preds, 1, quantile, probs = 0.975, na.rm = TRUE)


##Predicting Normal Growth Outcomes over 120-days ----
#load in weight-class groups
weight_groups <- readRDS("Inputs/p_24h_tbl.rds")
weight_groups <- weight_groups %>%
  mutate(Weight_g = Weight_kg *1000)  %>%
  select(Weight_g, Weight_Class) %>% 
  distinct(Weight_Class,.keep_all = TRUE)

glimpse(weight_groups)

# Function to calculate age from weight using inverted VBGF
# Weight = Winf * (1 - exp(-K * (Age - t0)))^3
# Solving for Age: Age = t0 - log(1 - (Weight/Winf)^(1/3)) / K
age_from_weight <- function(weight, Winf, K, t0) {
  t0 - log(1 - (weight / Winf)^(1/3)) / K
}

# Function to predict weight from age using VBGF
weight_from_age <- function(age, Winf, K, t0) {
  Winf * (1 - exp(-K * (age - t0)))^3
}

# Calculate initial age from initial weight
weight_groups <- weight_groups %>%
  mutate(
    initial_g = Weight_g,
    initial_age_years = age_from_weight(initial_g, Winf, K, t0),
    # Add 120 days (convert to years)
    final_age_years = initial_age_years + (120 / 365),
    # Predict final weight
    final_g = weight_from_age(final_age_years, Winf, K, t0),
    # Calculate growth
    growth_g = final_g - initial_g,
    percent_growth = (growth_g / initial_g) * 100
  )

# Display results
cat("\n=== 120-Day Growth Predictions ===\n")
print(weight_groups)

# Propagate bootstrap parameter uncertainty through 120-day predictions
boot_120 <- apply(bootTypical$coefboot, 1, function(p) {
  Winf_b <- p["Winf"]; K_b <- p["K"]; t0_b <- p["t0"]
  init_age_b  <- age_from_weight(weight_groups$initial_g, Winf_b, K_b, t0_b)
  final_age_b <- init_age_b + (120 / 365)
  weight_from_age(final_age_b, Winf_b, K_b, t0_b)
})
# boot_120 is a matrix: rows = weight classes, cols = bootstrap iterations
weight_groups$lower_final_g <- apply(boot_120, 1, quantile, probs = 0.025, na.rm = TRUE)
weight_groups$upper_final_g <- apply(boot_120, 1, quantile, probs = 0.975, na.rm = TRUE)

# Save weight_groups for use in BT FB4 model
saveRDS(weight_groups %>% select(Weight_Class, initial_g, final_g, lower_final_g, upper_final_g),
        "Inputs/weight_groups_vbgf.rds")
cat("\nSaved weight_groups to Inputs/weight_groups_vbgf.rds\n")

# Create growth curve plot
growth_plot <- ggplot() +
  # Bootstrap 95% CI band on curve
  geom_ribbon(data = pred_data, aes(x = Age, ymin = lower, ymax = upper),
              fill = "red", alpha = 0.2) +
  # Raw data points
  geom_point(data = growth_clean, aes(x = Age, y = RoundWeight),
             alpha = 0.4, size = 2, color = "steelblue") +
  # Fitted VBGF curve
  geom_line(data = pred_data, aes(x = Age, y = RoundWeight),
            color = "red", size = 1.2) +
  # Labels and theme
  labs(
       x = "Age (years)",
       y = "Weight (g)") +
  theme_classic() +
  theme(axis.line = element_line(linewidth = 1.5),
        axis.title = element_text(size = 18, face = "bold"),
        axis.text = element_text(size = 18, face = "bold"))

# Display the plot
print(growth_plot)
ggsave("Graphs/VBGF_growth_curve.png", growth_plot, width = 8, height = 6, dpi = 300)
