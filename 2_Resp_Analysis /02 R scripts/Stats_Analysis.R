######################################################
############################################################
#---------------------
  #STATISTICAL TESTS
#---------------------
  ############################################################

# Load packages
library(car)
library(pwr)
library(effectsize)
library(effsize)
library(WRS2)

#convert tempMeta to factor 
summary_df$TempMeta <- as.factor(summary_df$TempMeta)
summary_df$Treatment.x <- as.factor(summary_df$Treatment.x)

#check sample sizes. 
summary_df %>%
  group_by(TempMeta, Treatment.x) %>%
  summarise(n_fish = n_distinct(Fish.ID)) %>%
  arrange(TempMeta, Treatment.x)

# TempMeta Treatment.x n_fish
#  10       chase           12
# 10       chase + air     13
# 10       control         14
# 15       chase            8
# 15       chase + air      7
# 15       control          7

######################################################
############################################################
#---------------------
#BODY SIZE and MASS
# ---------------------
############################################################
str(data_SMR_Lab)

## Descriptive Summary Stats ###
summary_df$TL <- as.numeric(summary_df$TL)
summary_df$Weight.g <- as.numeric(summary_df$Weight.g)

#1) Body Size and TL 
# Overall body size and TL
All_size_stats <- summary_df %>%
  summarise(
    TL_mean = mean(TL, na.rm = TRUE),
    TL_sd = sd(TL, na.rm=TRUE),
    TL_range_min = min(TL, na.rm = TRUE),
    TL_range_max = max(TL, na.rm = TRUE),
    Weight_g_mean = mean(Weight.g, na.rm = TRUE),
    Weight_sd = sd(Weight.g, na.rm=TRUE),
    Weight_g_range_min = min(Weight.g, na.rm = TRUE),
    Weight_g_range_max = max(Weight.g, na.rm = TRUE)
  )
# Print the results
print(All_size_stats)


# Calculate summary statistics for each treatment group
Group_size_stats <- summary_df %>%
  group_by(Treatment.x, TempMeta) %>%   # Group by treatment
  summarise(
    TL_mean = mean(TL, na.rm = TRUE),
    TL_sd = sd(TL, na.rm=TRUE),
    TL_min = min(TL, na.rm = TRUE),
    TL_max = max(TL, na.rm = TRUE),
    Weight_g_mean = mean(Weight.g, na.rm = TRUE),
    Weight_sd = sd(Weight.g, na.rm=TRUE),
    Weight_g_min = min(Weight.g, na.rm = TRUE),
    Weight_g_max = max(Weight.g, na.rm = TRUE),
    
  )

# Print the results
print(Group_size_stats)

glimpse(summary_df)

## GRAPHING ##
# Create the boxplot
ggplot(summary_df, aes(x = Treatment.x, y = Weight.g, fill =(TempMeta))) +
  geom_boxplot(position = position_dodge(width = 0.8)) +
  labs(
    title = "Boxplot of Weight by Treatment and TempMeta",
    x = "Treatment",
    y = "Weight (g)",
    fill = "TempMeta"
  ) +
  theme_minimal()
######################################################
# INFERENTIAL STATS FOR TESTING TL AND WEIGHT
#make sure there are no differences in size between the groups
#we have unequal sample sizes and unequal variance across groups from levenes homogeneity test

# Homogeneity check
leveneTest(Weight.g ~ TempMeta * Treatment.x, data = summary_df)
#unequal variances here.

##WELCH's ANOVA
welch_anova <- t2way(Weight.g ~ TempMeta * Treatment.x, data = summary_df)
print(welch_anova)
#no signficant difference in weight between my groups

# Normality check
shapiro.test(residuals(m_weight))
qqnorm(residuals(m_weight)); qqline(residuals(m_weight))
#weight is normally distributed


##############################
######## Analyze Weight and SMR #######
####### ####### ####### ####### ####### ####### ####### ####### ####### ####### #######

#make sure correct order of treatment. 
summary_df$Treatment.x <- relevel(factor(summary_df$Treatment.x), ref = "control")

#1) MASS SPECIFIC MO2 and SMR/ MMR
# Fit the log-log linear model and test for treatment interaction
model <- lm(log10(MMR) ~ log10(Weight.g), data = summary_df)
summary(model)

par(mfrow = c(2, 2))  # Show 4 plots in one window
plot(model)
par(mfrow = c(1, 1))  # Reset to default

#plot by treatment
ggplot(summary_df, aes(x = log10(Weight.g), y = log10(MMR))) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", se = FALSE, color = "blue") +
  facet_wrap(~ Treatment.x) +
  labs(
    title = "Log-Log Plot of Mass Specific SMR vs Body Mass",
    x = "log10(Body Mass [g])",
    y = "log10(SMR [mg O2/kg/h])"
  ) 

#plot overall
ggplot(summary_df, aes(x = log10(Weight.g), y = log10(low10))) +
  geom_point() +  # Scatter plot
  geom_smooth(method = "lm", se = FALSE, color = "blue") +  # Line of best fit (linear model)
  labs(
    title = "Log Log Plot SMR vs Weight.g",  # Plot title
    x = "Weight (g)",  # X-axis label
    y = "SMR (mgO2/kg/h)"  # Y-axis label
  ) +
  theme_minimal()  # Optional: minimal theme for a clean look

######################################################
############################################################
#---------------------
#MO2 DATA
# ---------------------
############################################################
view(summary_df)

########### EPOC  ############

#1) EPOC Magnitude (ie: total EPOC)
epoc_mag_anova <- t2way(EPOC.magnitude ~ TempMeta * Treatment.x, data = summary_df)
print(epoc_mag_anova)

#Test shaprio-wilks test for normality 
shapiro.test(residuals(epoc_mag_anova)) 
qqnorm(residuals(epoc_mag_anova)); qqline(residuals(epoc_mag_anova))
#looks like assumptions of normality are not held!

# Levene's Test for homogeneity of variance
leveneTest(EPOC.magnitude ~ TempMeta *Treatment.x, data = summary_df)
#p>0.05. accept null. Data has homogeneous variance


#2) EPOC duration 
epoc_duration_aov <- aov(EPOC.duration ~ TempMeta * Treatment.x, data=summary_df)
summary(epoc_duration_aov)
#assumptions
leveneTest(EPOC.duration ~ TempMeta *Treatment.x, data = summary_df)
shapiro.test(residuals(epoc_duration_aov)) 
qqnorm(residuals(epoc_duration_aov)); qqline(residuals(epoc_duration_aov))

#3)Peak EPOC/SMR 
epoc_peak_aov <- aov(EPOC.peak ~ TempMeta * Treatment.x, data=summary_df)
summary(epoc_peak_aov)
#assumptions
leveneTest(EPOC.peak ~ TempMeta *Treatment.x, data = summary_df)
shapiro.test(residuals(epoc_peak_aov)) 
qqnorm(residuals(epoc_peak_aov)); qqline(residuals(epoc_peak_aov))
#significant difference in peak between the groups. 


######################################################################
########################################################
######### POWER ANALYSIS ##############
########################################################




