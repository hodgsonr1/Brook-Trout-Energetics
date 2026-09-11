######################################################
############################################################
#---------------------
#STATISTICAL TESTS
#---------------------

#Purpose: Statsitical tests for differences in treatments and temperature on 6 metrics of metabolic recovery 
#dependencies: summary_df output from Kenauk_Respirometry core script
############################################################

# Load packages
library(dplyr)
library(car) 
library(pwr)
library(effectsize)
library(effsize)
library(WRS2)
library(emmeans)
library(lme4)
library(lmerTest)
#--------------------------------------------------------
#load data 
summary_df <- readRDS("04_R_ouputs/summary_df.rds")

#create continous trial order 
summary_df <- summary_df %>%
  mutate(Trial = as.numeric(factor(paste(DateMeta, Trial), 
                                   levels = unique(paste(DateMeta, Trial)))))

#convert tempMeta to factor 
summary_df$TempMeta <- as.factor(summary_df$TempMeta)
summary_df$Treatment.x <- as.factor(summary_df$Treatment.x)
summary_df$TL <- as.numeric(summary_df$TL)
summary_df$Weight.g <- as.numeric(summary_df$Weight.g)


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
#-------------------------------------
## 1) Descriptive Stats 
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
print(All_size_stats)

#All treatment groups 
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
print(Group_size_stats)

## boxplot for body mass
ggplot(summary_df, aes(x = Treatment.x, y = Weight.g, fill = (TempMeta))) +
  geom_boxplot(position = position_dodge(width = 0.8)) +
  labs(
    title = "Boxplot of Weight by Treatment and TempMeta",
    x = "Treatment",
    y = "Weight (g)",
    fill = "TempMeta"
  ) +
  theme_minimal()
### test to make sure no difference in size across groups 
leveneTest(Weight.g ~ Treatment.x * TempMeta, data = summary_df)# unequal variances!
leveneTest(TL ~ Treatment.x * TempMeta, data=summary_df) 

#unequal variance in weight across treatments. 

#Weight ANOVA. 
options(contrasts = c("contr.sum", "contr.poly"))

lm_fit <- lm(Weight.g ~ TempMeta * Treatment.x, data = summary_df)
car::Anova(lm_fit, type = "III", white.adjust = "hc3")



#--------------------------
######################################################
############################################################
#---------------------
#MO2- EPOC Data 
# ---------------------
############################################################
#-------------------------------------
#Inferential stats
#-------------------------------------

#ANCOVAS
#--------------------

#2.1 magnitude of EPOC
mag_lm <- lmer(EPOC.magnitude ~ TempMeta * Treatment.x + Weight.g + (1| Trial), data=summary_df)
Anova(mag_lm, type = 3)
shapiro.test(residuals(mag_lm))
qqnorm(residuals(mag_lm)); qqline(residuals(mag_lm))
leveneTest(EPOC.magnitude ~ TempMeta * Treatment.x, data = summary_df)
summary(mag_lm)
#key results 
#---> no statistically significant findings!

#2.2 EPOC_duration
dur_lm <- lmer(EPOC.duration ~ TempMeta * Treatment.x + Weight.g + (1| Trial), data=summary_df)
Anova(dur_lm, type =3)
summary(dur_lm)
#check assumptions
shapiro.test(residuals(dur_lm)) 
qqnorm(residuals(dur_lm)); qqline(residuals(dur_lm))
leveneTest(EPOC.duration ~ TempMeta * Treatment.x, data = summary_df)

emmeans(dur_lm, pairwise ~ TempMeta, adjust = "tukey") 

contrasts(summary_df$TempMeta)

# #result --> longer recovery at 10C vs 15C for chase treatment 

#------------------------------------------------------------
#2.3 SMR
smr_lm <- lmer(low10 ~ TempMeta * Treatment.x + Weight.g + (1| Trial), data=summary_df)
Anova(smr_lm, type = 3)
summary(smr_lm)
contrasts()
shapiro.test(residuals(smr_lm)) #significant outliers here!
qqnorm(residuals(smr_lm));qqline(residuals(smr_lm))
leveneTest(low10 ~ TempMeta * Treatment.x, data=summary_df)
#key results 
#--> TempMeta p<0.001
#post hoc tests 
emmeans(smr_lm, pairwise ~ Treatment.x, adjust = "tukey")
#SMR is significantly larger at 15C compared to 10C. 
#------------------------------------------------------------
#2.4 peak EPOC

peak_lm <- lmer(EPOC.peak ~ TempMeta * Treatment.x + Weight.g + (1| Trial), data=summary_df)
Anova(peak_lm, type = 3) 
summary(peak_lm)
#check assumptions tests
shapiro.test(residuals(peak_lm)) 
qqnorm(residuals(peak_lm)); qqline(residuals(peak_lm))
leveneTest(EPOC.peak ~ TempMeta * Treatment.x, data = summary_df)
# Significant effect of temperature on the peak EPOC. 
#peak EPOC is significantly higher at 15C compared to 10C. 

#------------------------------------------------------------

#------------------------------------------------------------
#2.5 delta peak EPOC - SMR
delta_lm <- lmer(delta_peak ~ TempMeta * Treatment.x + Weight.g + (1| Trial), data=summary_df)
Anova(delta_lm, type=3)
summary(delta_lm)
shapiro.test(residuals(delta_lm))
qqnorm(residuals(delta_lm));qqline(residuals(delta_lm))
leveneTest(delta_peak ~ TempMeta * Treatment.x, data=summary_df)
#post hoc test
emmeans(delta_lm, pairwise ~ TempMeta, adjust = "tukey")
#key results 
#--> MO2 scope is higher at 15C compared to 10C.

#check order of temp (10 & 15) factors 
contrasts(summary_df$TempMeta)
contrasts(summary_df$Treatment.x)

#------------------------------------------------------------
#2.6 Ratio Peak EPOC/SMR
ratio_lm <- lmer(ratio ~ TempMeta * Treatment.x + Weight.g + (1| Trial) , data=summary_df)
Anova(ratio_lm, type=3)
summary(ratio_lm)
shapiro.test(residuals(ratio_lm))
qqnorm(residuals(ratio_lm));qqline(residuals(ratio_lm))
leveneTest(ratio ~ TempMeta * Treatment.x, data=summary_df)

emmeans(ratio_lm, pairwise ~ TempMeta, adjust = "tukey")

#key results 
#--> non significant





