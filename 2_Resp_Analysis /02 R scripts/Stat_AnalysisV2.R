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
library(emmeans)
#--------------------------------------------------------
#clean data 
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


#2.1 EPOC_duration
dur_lm <- lm(EPOC.duration ~ TempMeta * Treatment.x * Weight.g, data=summary_df)
Anova(dur_lm, type =3)
#check assumptions
shapiro.test(residuals(dur_lm)) 
qqnorm(residuals(dur_lm)); qqline(residuals(dur_lm))
leveneTest(EPOC.duration ~ TempMeta * Treatment.x, data = summary_df)
#Key Results
#--> Treatment - p=0.053* 
#--> Tempmeta:Treatment - p = 0.063*
#Effect of treatment on recovery duration (main effect)(Not interpretable because of interactions) 
#effect of temperature depends on treatment (interaction)
#effect of treatment depends on weight (interaction)

#POST-HOC - Tukey HSD.
#A) Do epoc durations at each treatment differ by temperature
emmeans(dur_lm, pairwise ~ TempMeta|Treatment.x, adjust = "tukey") 
#result --> longer recovery at 10C vs 15C for chase treatment 



#------------------------------------------------------------
#2.2 peak EPOC
peak_lm <- lm(EPOC.peak ~ TempMeta * Treatment.x * Weight.g, data=summary_df)
Anova(peak_lm, type = 3) 
#check assumptions tests
shapiro.test(residuals(peak_lm)) 
qqnorm(residuals(peak_lm)); qqline(residuals(peak_lm))
leveneTest(EPOC.peak ~ TempMeta * Treatment.x, data = summary_df)
#key results 
#--> TempMeta - p = 0.00617**
# Significant effect of temperature on the peak EPOC. 

#post hoc tests 
emmeans(peak_lm, pairwise ~ TempMeta, adjust = "tukey")
#peak EPOC is significantly higher at 15C compared to 10C. 

#------------------------------------------------------------
#2.3 magnitude of EPOC
mag_lm <- lm(EPOC.magnitude ~ TempMeta + Treatment.x + Weight.g, data=summary_df)
Anova(mag_lm, type = 3)
shapiro.test(residuals(mag_lm)) #significant but qq looks OK.
qqnorm(residuals(mag_lm)); qqline(residuals(mag_lm))
leveneTest(EPOC.magnitude ~ TempMeta * Treatment.x, data = summary_df)
#key results 
#---> no statistically significant findings!
#------------------------------------------------------------
#2.4 SMR
smr_lm <- lm(low10 ~ TempMeta * Treatment.x * Weight.g, data=summary_df)
Anova(smr_lm, type = 3)
shapiro.test(residuals(smr_lm)) #significant but qq looks OK
qqnorm(residuals(smr_lm));qqline(residuals(smr_lm))
leveneTest(low10 ~ TempMeta * Treatment.x, data=summary_df)
#key results 
#--> TempMeta p<0.001
#post hoc tests 
emmeans(smr_lm, pairwise ~ TempMeta, adjust = "tukey")
#SMR is significantly larger at 15C compared to 10C. 

#------------------------------------------------------------
#2.5 delta peak EPOC - SMR
delta_lm <- lm(delta_peak ~ TempMeta + Treatment.x + Weight.g, data=summary_df)
Anova(delta_lm, type=3)
shapiro.test(residuals(delta_lm))
qqnorm(residuals(delta_lm));qqline(residuals(delta_lm))
leveneTest(delta_peak ~ TempMeta * Treatment.x, data=summary_df)
#post hoc test
emmeans(delta_lm, pairwise ~ TempMeta, adjust = "tukey")
#key results 
#--> MO2 scope is higher at 15C compared to 10C.



#------------------------------------------------------------
#2.6 Ratio Peak EPOC/SMR
ratio_lm <- lm(ratio ~ TempMeta * Treatment.x * Weight.g, data=summary_df)
Anova(ratio_lm, type=3)
shapiro.test(residuals(ratio_lm))
qqnorm(residuals(ratio_lm));qqline(residuals(ratio_lm))
leveneTest(ratio ~ TempMeta * Treatment.x, data=summary_df)
#key results 
#--> non significant





