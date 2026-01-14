### Title: 
#MNR Growth Data Extraction ###

##Description: 
#load in brook trout age-growth data from MNR long term monitoring program. 
#use to build growth over time relationship using VGBF

#Ryan Hodgson 
#Nov 29 2025

getwd()

#load packages
library(dplyr)
library(FSA)
library(ggplot2)

#load data
growth <- read.csv("Inputs/BsM_BTdata.csv")

head(growth)

#clean data ----
# Select relevant columns and remove missing values
growth_clean <- growth %>%
  select(Age, RoundWeight) %>%
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


#fit VBGF model ----
# Examine max values to inform starting parameters
max_age <- max(growth_clean$Age)
max_weight <- max(growth_clean$RoundWeight)
cat("\nMax observed age:", max_age, "years\n")
cat("Max observed weight:", max_weight, "g\n")

# Set manual starting values based on biological knowledge
# Winf should be larger than max observed weight (asymptotic value)
# K typically ranges from 0.1-0.5 for fish
# t0 typically slightly negative
sv <- list(Winf = max_weight * 1.2,  # Set Winf 20% higher than max observed
           K = 0.3,                    # Moderate growth rate
           t0 = -0.5)                  # Typical negative t0

cat("\nManual starting parameter values:\n")
cat("Winf:", sv$Winf, "g\n")
cat("K:", sv$K, "\n")
cat("t0:", sv$t0, "\n")

# Fit the VBGF model using nonlinear least squares
fit <- nls(RoundWeight ~ Winf * (1 - exp(-K * (Age - t0)))^3,
           data = growth_clean,
           start = sv)

# Display model summary
cat("\n=== VBGF Model Summary ===\n")
summary(fit)

# Extract parameters
params <- coef(fit)
Winf <- params["Winf"]
K <- params["K"]
t0 <- params["t0"]

cat("\n=== Fitted Parameters ===\n")
cat("Winf (asymptotic weight):", round(Winf, 2), "g\n")
cat("K (growth coefficient):", round(K, 4), "\n")
cat("t0 (theoretical age at weight=0):", round(t0, 2), "years\n")

#model diagnostics ----
# Calculate residuals and fitted values
growth_clean$fitted <- predict(fit)
growth_clean$residuals <- residuals(fit)

# Calculate goodness-of-fit metrics
SS_res <- sum(growth_clean$residuals^2)
SS_tot <- sum((growth_clean$RoundWeight - mean(growth_clean$RoundWeight))^2)
R_squared <- 1 - (SS_res / SS_tot)
RMSE <- sqrt(mean(growth_clean$residuals^2))

cat("\n=== Model Fit Statistics ===\n")
cat("R-squared:", round(R_squared, 4), "\n")
cat("RMSE:", round(RMSE, 2), "g\n")
cat("Mean absolute error:", round(mean(abs(growth_clean$residuals)), 2), "g\n")

# Create diagnostic plots
# 1. Residuals vs Fitted
p1 <- ggplot(growth_clean, aes(x = fitted, y = residuals)) +
  geom_point(alpha = 0.5, color = "steelblue") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  geom_smooth(se = FALSE, color = "black", linewidth = 0.5) +
  labs(title = "Residuals vs Fitted Values",
       x = "Fitted Weight (g)",
       y = "Residuals (g)") +
  theme_bw()

# 2. Residuals vs Age
p2 <- ggplot(growth_clean, aes(x = Age, y = residuals)) +
  geom_point(alpha = 0.5, color = "steelblue") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  geom_smooth(se = FALSE, color = "black", linewidth = 0.5) +
  labs(title = "Residuals vs Age",
       x = "Age (years)",
       y = "Residuals (g)") +
  theme_bw()

# 3. Histogram of residuals
p3 <- ggplot(growth_clean, aes(x = residuals)) +
  geom_histogram(bins = 30, fill = "steelblue", color = "black", alpha = 0.7) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  labs(title = "Distribution of Residuals",
       x = "Residuals (g)",
       y = "Frequency") +
  theme_bw()

# Display diagnostic plots
cat("\nDisplaying diagnostic plots...\n")
print(p1)
print(p2)
print(p3)

#VGBF visualization ----
# Create prediction data for smooth curve
age_pred <- seq(min(growth_clean$Age), max(growth_clean$Age), length.out = 100)
weight_pred <- predict(fit, newdata = data.frame(Age = age_pred))

pred_data <- data.frame(Age = age_pred, RoundWeight = weight_pred)

# Create growth curve plot
growth_plot <- ggplot() +
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

# Save weight_groups for use in BT FB4 model
saveRDS(weight_groups %>% select(Weight_Class, initial_g, final_g),
        "Inputs/weight_groups_vbgf.rds")
cat("\nSaved weight_groups to Inputs/weight_groups_vbgf.rds\n")
