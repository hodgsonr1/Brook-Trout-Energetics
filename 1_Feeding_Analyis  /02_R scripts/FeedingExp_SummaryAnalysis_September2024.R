########### FEEDING EXPERIMENT 2024 #############
###### DESCRIPTIVE STATS AND PRELIMENARY SUMMARY ###########################

# load in the data
Data <- read.csv("FeedingExpTrialData_Kenauk2024V2.csv")
#load packages
library(dplyr)
library(ggplot2)
# summarize and count sample size for each of the treatments.
table(Data$Treatment)

## CLEAN DATA ##
# Remove rows with NA in the Feed.Delivery.Number column
Data_clean <- Data[!is.na(Data$Feed.Delivery.Number), ]

####### CALCULATING MEAN/Median FEEDING LATENCY FOR EACH TREATMENT ################
############################################################################################

str(Data)
numeric(Data$Weight..g.)
mean(Data$WTotal.Length..mm.)

# Calculate the summary statistics
Data_summary1 <- Data_clean %>%
  group_by(Treatment) %>%
  summarize(
    mean_value = mean(Feed.Delivery.Number, na.rm = TRUE),
    median_value = median(Feed.Delivery.Number, na.rm = TRUE),
    sd = sd(Feed.Delivery.Number, na.rm=TRUE),
    min = min(Feed.Delivery.Number, na.rm=TRUE),
    max = max(Feed.Delivery.Number, na.rm=TRUE)
    
  )

print(Data_summary1)

######### calculate standard error (SEM) and Deviation (SD) stats #
Error_summary_data <- Data_clean%>%
  group_by(Treatment == "control") %>%
  summarize(
    mean = mean(Feed.Delivery.Number, na.rm = TRUE),
    sd = sd(Feed.Delivery.Number, na.rm = TRUE),
    n = n(),
    se = sd / sqrt(n),  # Standard Error of the Mean
    .groups = 'drop'
  )
print(Error_summary_data)
############################## ############################## ##############################
################ GRAPHING OUT FEEDING LATENCY DATA ##############################
############################## ############################## ##############################

# create a boxplot of feeding latencies 

ggplot(Data_clean, aes(x = Treatment, y = Feed.Delivery.Number)) + ## note used feed.delivery.number because it is equivalent to feeding latency in hours
  geom_boxplot() +
  stat_boxplot(geom = "errorbar", 
               width = 0.25) +
  labs(title = "Feeding Latency by Treatment", 
       x = "Treatment", 
       y = "Feeding Latency (hrs)") +
  scale_x_discrete(labels = c("chase" = "chase", 
                              "chase_air" = "chase & air", 
                              "control_A" = "control (A)",
                              "control_B" = "control (B)" )) +  # Change x-axis labels
  theme_minimal() +
  theme (panel.grid.major.y = element_blank(),  # Remove horizontal major gridlines
                panel.grid.minor.y = element_blank(),  # Remove horizontal minor gridlines
                axis.line = element_line(),
         panel.grid.major.x = element_blank(),  # Optional: Remove vertical gridlines
         panel.grid.minor.x = element_blank(),  # Optional: Remove vertical minor gridlines,
         plot.title = element_text(hjust = 0.5, size =16),
         axis.title.x = element_text(size=15),
         axis.title.y = element_text(size=15),
         axis.text.x = element_text(size = 12, color = "black"),  # Set font size and color for x-axis labels
         axis.text.y = element_text(size = 12, color = "black"))  # Keep axis lines)

###################################################################################################
############################## GRAPHING FEEDING LATENCY BY Mass and TL ###########################
######################################################################################################

### Total Mass
ggplot(Data_clean, aes(x= Weight..g., y=Feed.Delivery.Number, color = Treatment))+
       geom_point( size=3)+
         labs(title = "Distribution of Feeding Latency by Weight ", x = " Weight (g)", y = "Feeding latency (hrs)")+
  theme_minimal() +
  theme (panel.grid.major.y = element_blank(),  # Remove horizontal major gridlines
         panel.grid.minor.y = element_blank(),  # Remove horizontal minor gridlines
         axis.line = element_line(),
         panel.grid.major.x = element_blank(),  # Optional: Remove vertical gridlines
         panel.grid.minor.x = element_blank(),  # Optional: Remove vertical minor gridlines,
         plot.title = element_text(hjust = 0.5, size =16),
         axis.title.x = element_text(size=15),
         axis.title.y = element_text(size=15),
         axis.text.x = element_text(size = 12, color = "black"),  # Set font size and color for x-axis labels
         axis.text.y = element_text(size = 12, color = "black"))  # Keep axis lines)
       
# Total Length
ggplot(Data_clean, aes(x= Total.Length..mm., y=Feed.Delivery.Number, color = Treatment))+
  geom_point( size=3)+
  labs(title = "Distribution of Feeding Latency by Total Length", x = " Total Length (mm)", y = "Feeding latency (hrs)")+
  theme_minimal() +
  theme (panel.grid.major.y = element_blank(),  # Remove horizontal major gridlines
         panel.grid.minor.y = element_blank(),  # Remove horizontal minor gridlines
         axis.line = element_line(),
         panel.grid.major.x = element_blank(),  # Optional: Remove vertical gridlines
         panel.grid.minor.x = element_blank(),  # Optional: Remove vertical minor gridlines,
         plot.title = element_text(hjust = 0.5, size =16),
         axis.title.x = element_text(size=15),
         axis.title.y = element_text(size=15),
         axis.text.x = element_text(size = 12, color = "black"),  # Set font size and color for x-axis labels
         axis.text.y = element_text(size = 12, color = "black"))  # Keep axis lines)

