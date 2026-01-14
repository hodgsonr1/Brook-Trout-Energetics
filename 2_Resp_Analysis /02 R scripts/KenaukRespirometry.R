##################################################################
# Title: Respirometry Analysis Kenauk 2024/25

# Purpose:This script transforms and combines respriometry trial data with meta data. Several calculations are performed including MO2, EPOC, SMR and MMR. 
#         Data used in from experients at Kenauk in 2024 and 2025, examining the metabolic costs of simulated catch and release angling events, on hatchery brook trout.  

# Date created: October 18 2024
# Date last modified: Nov 7 2025 ------------------------------
RStudio.Version()

#
#Author: Ryan Hodgson
#Contributions from: Paul Bznonek, Christian Bihun and Jake Brownscombe
##################################################################


#Load packages
library(tidyverse) #ggplot, dplyr, etc.
library(plotly) #interactive plots
library(rio) #Import lists of files
library(readxl) #Real excel files
library(mclust) #gaussian mixture model (GMM) for type of SMR calc from Chabott (2016)
library(ggsignif) #for adding significance brackets to plots
library(quantreg) #used for smoothing epoc curves using quantile smoothing spleen
library(SparseM) # used for creating EPOC curves and smoothing. 
library(Hmisc) #again used for creating EPOC curves and smoothing. 


#Set global environment
theme_set(theme_classic()) #Set global ggplot theme
options(scipen=999) #Remove scientific notation

#check working directory 
getwd()


#####Read in the Data#########################################----
#-------------------------------------------------------------# 
###Get the Meta data----#
Lab_MetaData <- read_xlsx("01 Good Raw Data/MetaData_BrookTrout_KenaukResp_Combined.xlsx") 
names(Lab_MetaData) <- c("Date", "Trial", "Fish.ID","BackgroundVsTrial", "Holding Tank", "Temp", "Treatment", "Chamber", "Channel",                
                                  "Trial Start", "Trial End", "Flush Timing", "TL", "Weight.g", "Notes", "source")
Lab_MetaData <- Lab_MetaData %>% 
  mutate(channel = paste0("Ch", Chamber), #create new column that concatenates Ch and chamber # (ie: Ch1 etc)
         Date=as.Date(Date, tryFormats = c("%d-%m-%Y")), #reformats date values
         nrow=nrow(.), #adds n columns with values of 99 (# of lines in dataframe)
         Trial_Start = as.POSIXct(`Trial Start`), #adds new row with start and ends using POSIXCt format. ie: seconds formatting.
         Trial_End = as.POSIXct(`Trial End`), #
  )
str(Lab_MetaData)

### Load in Trial Data ----- #
#Make the clean list of raw files
files <- append(x = c(), #Make a new vector to build from
                list.files(path = "01 Good Raw Data", #Look in location...
                           pattern = ".csv", recursive=T)) #Pull out `.csv` files from folders
files <- grep("raw", files, value=TRUE) #Remove files that don't contain "raw"
files <- paste0("01 Good Raw Data/", files) #return files names to file path

#Import the data
#NOTE: Column order must be maintained across files
data_SMR_Lab_Raw <- import_list(files, rbind = TRUE, rbind_label = "source", #Make a single table
                                skip=1)  #Skip header to avoid syntax variation
rm(files)

print(unique(data_SMR_Lab_Raw$source))
#remove irrelevant columns
data_SMR_Lab_Raw <- data_SMR_Lab_Raw[,c(1,2,3,6,7,10,13,16,48)] #only keep these columns.
names(data_SMR_Lab_Raw) <- c("datetime", "Phase", "Barometric", "Temp", 
                             "Ch1", "Ch2", "Ch3", "Ch4", "source")


#reformat data and determine closed periods
data_SMR_Lab_Raw <- data_SMR_Lab_Raw %>% 
  pivot_longer(cols=c("Ch1":"Ch4"), names_to="channel", values_to="O2") %>% #lengthen data
  mutate(datetime = as.POSIXct(datetime, tz="EST", format="%Y-%m-%d/%I:%M:%S %p"), #Date formatting
         Date = as.Date(datetime),# add new row with just date (not time)
         period = as.numeric(substr(Phase, 2,4))-1, #Grab loop number from phase column. 
         type2 = grepl("M", Phase), #Identify closed periods (return TRUE FALSE list if M appears in phase)
         type2 = case_when(type2==TRUE ~ "Close", #change TRUE to close otherwise flush.
                           TRUE ~ "Flush"),
  ) %>%
  mutate(across(where(is.character), ~as.factor(.)))

#Only keep closed loops
data_SMR_Lab <- data_SMR_Lab_Raw %>% 
  filter(period > 0, 
         !is.nan(O2),
         type2 == "Close"
  ) #Drop first loop

##############

## Add metadata, noting that the days wont match-up for overnight trials
# ---------------------------_ #######

#combine the two dataframes 
data_SMR_Lab <- left_join(data_SMR_Lab, Lab_MetaData, 
                          by = c("source","channel"), relationship = "many-to-many")


# remove all whitespaces
data_SMR_Lab <- data_SMR_Lab %>%
  mutate(datetime = str_trim(datetime),
         Trial_Start = str_trim(Trial_Start),
         Trial_End = str_trim(Trial_End))

# keep only selected columns and rename 
data_SMR_Lab <-data_SMR_Lab %>%
  select(Fish.ID, Trial, TempMeta=Temp.y, DateMeta=Date.y, DateActual=Date.x, datetime,         
         Treatment, BackgroundVsTrial, Chamber, Channel,              
         Trial_Start, Trial_End, Phase, TL, Weight.g,   
         TempActual=Temp.x, O2, period, type2, Barometric, Notes, source)



##################################################
### filtering the combined data 
#___________________###
# filter the new file to keep only values within start and end times
data_SMR_Lab <-data_SMR_Lab %>%
  filter(datetime >= Trial_Start,
         datetime <= Trial_End) 
#check to see if all files remain
print(unique(data_SMR_Lab$source))
gc() #Free unused memory

#########################################################

#ensure correct formatting of date time
data_SMR_Lab <- data_SMR_Lab %>%
  mutate(datetime = as.POSIXct(datetime, format = "%Y-%m-%d %H:%M:%S"),
         Trial_Start = as.POSIXct(Trial_Start, format = "%Y-%m-%d %H:%M:%S")
  )


#Trim start and end times of each loop by 30s
df_timebuffer <- data_SMR_Lab %>% 
  group_by(Fish.ID, period) %>% 
  summarise(Trial.Start = min(datetime),
            Trial.End = max(datetime),
            Trial.Start.buffer = Trial.Start + 30,
            Trial.End.buffer = Trial.End - 30,
            Reference.Temp = first(TempActual),#calculate reference temp
            ) %>%   
  #calculate trial duration
  group_by(Fish.ID) %>% # group by Fish.ID to calculate durations
  # Order by Trial.Start to calculate duration in the correct order
  arrange(Fish.ID, Trial.Start) %>% 
  # Create a new column with the time difference between consecutive periods
  mutate(Duration = difftime(Trial.Start, lag(Trial.Start), units = "hours")) %>% 
  # Replace NA in the first period with 0 (since there's no previous period for the first one)
  mutate(Duration = ifelse(is.na(Duration), 0, Duration)) %>%
  # Calculate cumulative sum of the durations
  mutate(Cumulative.Duration = cumsum(Duration)) %>%
  ungroup()  # Ungroup after calculations

#add data back in
data_SMR_Lab <- left_join(data_SMR_Lab, df_timebuffer, 
                          by = c("Fish.ID", "period")) %>% 
  filter(datetime>=Trial.Start.buffer & datetime<=Trial.End.buffer) 

#Convert Fish.ID to numeric for modeling
data_SMR_Lab <- data_SMR_Lab %>% 
  mutate(Fish.ID.n = as.numeric(as.factor(Fish.ID)),
         datetime.n = as.numeric(datetime))
#same with time
data_SMR_Lab <- data_SMR_Lab %>%
group_by(Fish.ID, period) %>%
  mutate(datetime.n = row_number()) %>%
  ungroup()



#####################################################################
#----------- TEMPERATURE ADJUSTMENT CODE  -----------##########
####
## Make function to calculate full saturation
oxygen_saturation <- function(T) {
  14.621 - 0.41022 * T + 0.007991 * T^2 - 0.000077774 * T^3
}
#make function to calculate DO%
DO_percentage <- function(measured_DO, T) {
  saturation_DO <- oxygen_saturation(T)
  DO_percent <- (measured_DO / saturation_DO) * 100
  return(DO_percent)
}
#make sure data are numeric
data_SMR_Lab$Reference.Temp <-
  as.numeric(as.character(data_SMR_Lab$Reference.Temp))
data_SMR_Lab$TempActual <-
  as.numeric(as.character(data_SMR_Lab$TempActual))
data_SMR_Lab$O2 <-
  as.numeric(as.character(data_SMR_Lab$O2))

#temperature correction 
data_SMR_Lab$full_sat_reference <- oxygen_saturation(data_SMR_Lab$Reference.Temp) # Apply function to reference temperature 
data_SMR_Lab$full_sat_actual <- oxygen_saturation(data_SMR_Lab$TempActual) # And to actual temperatures
data_SMR_Lab$sat_diff<-data_SMR_Lab$full_sat_reference -  data_SMR_Lab$full_sat_actual#Calculate the difference in O2 saturation between the reference and each timepoint
data_SMR_Lab$O2_adjusted <- data_SMR_Lab$O2 - data_SMR_Lab$sat_diff  #Adjust O2

#
# Calculate DO_percent_actual based on actual O2 and temperature
data_SMR_Lab$DO_percent_actual <- DO_percentage(data_SMR_Lab$O2, data_SMR_Lab$TempActual) 
data_SMR_Lab$DO_percent_adjusted <- DO_percentage(data_SMR_Lab$O2_adjusted, data_SMR_Lab$TempActual)
data_SMR_Lab$DO_percent_reference_check <- DO_percentage(data_SMR_Lab$full_sat_reference, data_SMR_Lab$Reference.Temp)

#ensure that function is calculating DO% and Mg/L correctly. ie: they match up.
all_close_to_100 <- all(abs(data_SMR_Lab$DO_percent_reference_check - 100) < 0.001)
if (all_close_to_100) {
  print("All reference DO% values are correctly calculated as 100% saturation.")
} else {
  print("Warning: Some reference DO% values are not exactly 100%!")
}


#####BACKGROUND TRIALS ###################################----
#-------------------------------------------------------------# 
str(data_SMR_Lab)

#make data frame for background trials
data_background <- data_SMR_Lab %>% 
  filter(BackgroundVsTrial == "Background")
               
#filter the background trials to 1hr and update periods to calculate slopes every 10 mins #
data_background <- data_background %>%
  group_by(Fish.ID) %>%
  filter(datetime >= Trial_Start & datetime <= Trial_Start + 60 * 60) %>%  # Filter for 1 hour after Trial_Start
  mutate(
    time_diff_secs = as.numeric(difftime(datetime, Trial_Start, units = "secs")),  # Calculate time difference
    period = floor(time_diff_secs / 600) + 1  # Update period based on 600s intervals
  ) %>%
  ungroup()

#filter to keep only manually selected trials which look reasonable (ie: no upward spikes, and irregularities in O2)
data_background <- data_background %>%
  filter(
    Fish.ID %in% c("2024-08-04- BLANK 08", #trial 1 10c
                   "2024-08-11 - BLANK 01" ,#trial 4 10c
                   "2024-08-11 - BLANK 02",
                   "2024-08-11 - BLANK 03",
                   "2024-08-11 - BLANK 04", 
                   "2024-08-17 - BLANK 01",#trial 7 10c
                   "2024-08-17 - BLANK 02",
                  "2024-08-17 - BLANK 03",
                   "2024-08-17 - BLANK 04",#trial 8 10c
                  "2024-08-19 - BLANK 01",
                  "2024-08-19 - BLANK 02",
                  "2024-08-19 - BLANK 03",
                  "2024-08-19 - BLANK 04",
                  "2024-08-22 - BLANK 05",#trial 9 10c
                  "2024-08-22 - BLANK 06",
                  "2024-08-22 - BLANK 07",
                  "2024-08-22 - BLANK 08",
                  "2024-08-25 - BLANK 04", #trial 11 10c
                  "2024-08-25 - BLANK 05",
                  "2024-08-25 - BLANK 06",
                  "2025-06-02-BLANK 01", #trial 1 15c
                  "2025-06-02-BLANK 02",
                  "2025-06-02-BLANK 03",
                  "2025-06-07-BLANK 04",#trial 2 15c
                  "2025-06-07-BLANK 05",
                  "2025-06-07-BLANK 06", #trial 5 15c
                  "2025-06-11-BLANK 01",
                  "2025-06-11-BLANK 02",
                  "2025-06-11-BLANK 03"
                   ) 
  )

### ANALYZE AND CALCULATE SLOPES

#now calculate slopes
#Loop through backgrounds and get summary lm results
loopB <- c()

# Loop through each Fish.ID
for (i in min(data_background$Fish.ID.n):max(data_background$Fish.ID.n)) {
  loopA <- filter(data_background, Fish.ID.n == i)  # Filter by Fish.ID.n
  for (j in unique(loopA$period)) {  
    tryCatch({
      loopA_period <- filter(loopA, period == j)  # Filter by period 
      
      # Linear model on O2 adjusted by datetime.n
      lm <- lm(O2_adjusted ~ datetime.n, data = loopA_period)
      
      # Collect the outputs
      outputs <- cbind(
        Trial = loopA_period$Trial[1],
        TempMeta = loopA_period$TempMeta[1],
        Fish.ID = loopA_period$Fish.ID[1],
        O2max = max(loopA_period$O2_adjusted),
        O2min = min(loopA_period$O2_adjusted),
        O2mean = mean(loopA_period$O2_adjusted),
        O2sd = sd(loopA_period$O2_adjusted),
        O2intercept = lm$coefficients['(Intercept)'],
        O2slope = -(lm$coefficients['datetime.n']),
        r2 = summary(lm)$r.squared
      )
      
      # Append outputs to loopB
      loopB <- rbind(loopB, outputs)
    }, error = function(e) { NULL })  # Handle any errors gracefully
  }
}
# Convert loopB to a dataframe
data_background_results <- as.data.frame(loopB); rm(loopB)
row.names(data_background_results) <- NULL

# Reformat the results
data_background_results <- data_background_results %>% 
  mutate(
    O2intercept = as.numeric(as.character(O2intercept)),
    O2slope = as.numeric(as.character(O2slope)),
    O2max = as.numeric(as.character(O2max)),
    O2min = as.numeric(as.character(O2min)),
    O2mean = as.numeric(as.character(O2mean)),
    O2sd = as.numeric(as.character(O2sd)),
    r2 = as.numeric(as.character(r2)) 
  )  %>%
  filter(( O2slope > 0))#filter out positive slopes


#get average stats for each temp to later apply... 

background_trial_avg <- data_background_results %>%
  dplyr::group_by(TempMeta) %>%  # Group by Trial
  dplyr::summarize(
    B_O2max_avg = mean(O2max, na.rm = TRUE),
    B_O2min_avg = mean(O2min, na.rm = TRUE),
    B_O2mean_avg = mean(O2mean, na.rm = TRUE),
    B_O2sd_avg = mean(O2sd, na.rm = TRUE),
    B_O2intercept_avg = mean(O2intercept, na.rm = TRUE),
    B_O2slope_avg = mean(O2slope, na.rm = TRUE),
    B_r2_avg = mean(r2, na.rm = TRUE)
  )

#################################################
##### TRIAL DATA######
#------------------------------
############

#make a numeric value column for linear model fitting for each fishID and period:
data_SMR_Lab <- data_SMR_Lab %>%
  group_by(Fish.ID, period) %>%
  mutate(datetime.n = row_number()) %>%
  ungroup()

#Loop through experiments and get summary lm results
#this loop goes through each fish ID and period, selects the data and fits a linear model to it, pulling out a range of O2 values,
#as well as the intercept, slope, r2 of the model, which are the main things you need to calculate O2 consumption
loopC <- c()
for (i in min(data_SMR_Lab$Fish.ID.n, na.rm = T):max(data_SMR_Lab$Fish.ID.n, na.rm = T)){
  loopA <- filter(data_SMR_Lab, Fish.ID.n==i)
  for (j in unique(loopA$period)){
    tryCatch({
      loopB <- filter(loopA, period==j)
      lm <- lm(O2_adjusted~datetime.n, data=loopB)
      outputs <- cbind(Fish.ID=loopB$Fish.ID[1],
                       Treatment=loopB$Treatment[1],
                       Trial = loopB$Trial[1],
                       period=loopB$period[1],
                       O2max=max(loopB$O2_adjusted),
                       O2min=min(loopB$O2_adjusted),
                       O2mean=mean(loopB$O2_adjusted),
                       O2sd=sd(loopB$O2_adjusted),
                       O2intercept=lm$coefficients['(Intercept)'], 
                       O2slope=-(lm$coefficients['datetime.n']),
                       r2=summary(lm)$r.squared)
      loopC <- as.data.frame(rbind(loopC, outputs))
    }, error=function(e){NULL})
  }}; rm(loopA); rm(loopB); rm(lm); rm(outputs); rm(i); rm(j)

head(loopC)

#Paste and reformat results
data_SMR_Lab_results <- as.data.frame(loopC); rm(loopC)
row.names(data_SMR_Lab_results) <- NULL
data_SMR_Lab_results <- data_SMR_Lab_results %>% 
  mutate(O2intercept = as.numeric(as.character(O2intercept)),
         O2slope = as.numeric(as.character(O2slope)),
         O2max = as.numeric(as.character(O2min)),
         O2min = as.numeric(as.character(O2min)),
         O2mean = as.numeric(as.character(O2mean)),
         O2sd = as.numeric(as.character(O2sd)),
         r2 = as.numeric(as.character(r2)),
         period = as.numeric(as.character(period)))


print(unique(data_SMR_Lab_results$Fish.ID))
#----------------------------#
###Convert slopes to O2/L/h
#----------------------------#

str(data_SMR_Lab)
str(data_SMR_Lab_results)
str(background_trial_avg)
background_trial_avg$TempMeta <- as.numeric(background_trial_avg$TempMeta)

#----- Calculate Mass specific metabolism
data_SMR_Lab_metabolism <- 
  left_join(data_SMR_Lab_results, #Merge results and metadata
            unique(select(data_SMR_Lab,Treatment, Fish.ID, DateMeta, TempMeta, Chamber,  TL, Weight.g)),
            by=c("Fish.ID")) %>% 
  left_join(., unique(select(data_SMR_Lab, Fish.ID, period, Trial.Start, Cumulative.Duration)), by=c("Fish.ID", "period")) %>% 
  # Join with background_trial_avg to apply background corrections
  left_join(., background_trial_avg, by = "TempMeta") %>%
  mutate(Weight.g = as.numeric(Weight.g),
         Volume = case_when(Chamber == "1" ~ (6.615), 
                            Chamber == "2" ~ (6.650),
                            Chamber == "3" ~ (6.570),
                            Chamber == "4" ~ (6.690)),
         #Metabolism equation
         O2.L.h_Raw = O2slope*3600*(Volume-(Weight.g/1000)),
         O2.L.h_Background = B_O2slope_avg*3600*Volume,
         O2.L.h_Corrected = O2.L.h_Raw - O2.L.h_Background,
         O2.Kg.L.h = O2.L.h_Corrected/(Weight.g/1000),
         O2.Kg.L.h_raw = O2.L.h_Raw/(Weight.g/1000),
         BackgroundO2Ratio = O2.L.h_Background/O2.L.h_Raw,
         
         #Identify 'valid' resp estimates
         SMR.valid = ifelse(test = (O2slope > 0 & #Slope is positive (O2 declines over time)
                                      r2 > 0.9   &  #R2 cutoff of 0.9
                                      O2.Kg.L.h > 0 & #Oxygen consumption cant be negative
                                      !Fish.ID %in% grep("Background|Blank", Fish.ID, value=TRUE)), #ignore background resp
                            yes=TRUE,
                            no=FALSE),
  ) %>% filter(!is.na(O2slope))



#----------------------------#
###Clean up the Data
#----------------------------#

#remove background rows

# Remove rows where SMR.valid is FALSE. ie: where R^2 <0.9
data_SMR_Lab_metabolism <- data_SMR_Lab_metabolism %>%
  filter(SMR.valid == TRUE)

#filter data down to 24hrs only 
data_SMR_Lab_metabolism <- data_SMR_Lab_metabolism[data_SMR_Lab_metabolism$Cumulative.Duration <= 24.00, ]

#remove 2025-06-10-003- looks like chamber was leaking. Only 3 measurement with R^2 > 0.9
data_SMR_Lab_metabolism <- data_SMR_Lab_metabolism %>%
  filter(Fish.ID != "2025-06-10-003")

#----------------------------#
###Calculating SMR
#----------------------------#
########################################################################################
#------------------------------------------


#1) D.Chabott 2016 SMR method comparisons 
  #--------------------------#

  # Initialize the dataframe to store the results
SMR_results <- data.frame(
    Fish.ID = numeric(0),
    mlnd = numeric(0),
    CVmlnd = numeric(0),
    low10 = numeric(0),
    low10pc = numeric(0),
    quant_10 = numeric(0),  # Quantile for 10% (q = 0.1)
    quant_15 = numeric(0),  # Quantile for 15% (q = 0.15)
    quant_20 = numeric(0),  # Quantile for 20% (q = 0.2)
    quant_25 = numeric(0),  # Quantile for 25% (q = 0.25)
    quant_30 = numeric(0)   # Quantile for 30% (q = 0.3)
)


# Loop through each unique Fish.ID
for (fish in unique(data_SMR_Lab_metabolism$Fish.ID)) {
  
  # Subset data for the current fish
  fish_data <- subset(data_SMR_Lab_metabolism, Fish.ID == fish)
  
  # Filter for valid data (SMR.valid == TRUE)
  valid_data <- subset(fish_data, SMR.valid == TRUE)
  
  # Proceed if there is valid data for the fish
  if (nrow(valid_data) > 1) {
    
    # Get the MO2 data for the current fish
    Y <- valid_data$O2.Kg.L.h
    
    # Fit a Gaussian Mixture Model (GMM) to the MO2 data (1 to 4 components)
    the.Mclust <- Mclust(Y, G = 1:4)
    
    # Get the classification of each data point
    cl <- the.Mclust$classification
    
    # Identify the cluster that likely represents the SMR (using the provided logic)
    cl2 <- as.data.frame(table(cl))
    cl2$cl <- as.numeric(levels(cl2$cl))
    valid <- cl2$Freq >= 0.1 * length(Y)  # Cluster is valid if it contains at least 10% of the data
    the.cl <- min(cl2$cl[valid])
    
    # Get the MO2 values from the identified SMR cluster
    left.distr <- Y[the.Mclust$classification == the.cl]
    
    # Calculate the mean of the lowest normal distribution (MLND)
    mlnd <- the.Mclust$parameters$mean[the.cl]
    
    # Calculate the Coefficient of Variation (CV) for the SMR cluster
    CVmlnd <- sd(left.distr) / mlnd * 100
    
    # Calculate quantiles for the data (using the specified q values)
    q <- c(0.1, 0.15, 0.2, 0.25, 0.3)
    quant <- quantile(Y, q)
    
    # Calculate the average of the lowest 10 values
    u <- sort(Y)
    low10 <- mean(u[1:10])
    
    # Calculate the low 10% average, excluding outliers (method from Herrmann & Enders 2000)
    low10pc <- mean(u[6:(5 + round(0.1 * (length(u) - 5)))])  # Exclude first 5 outliers
    
    # Store the results for this fish
    SMR_results <- rbind(SMR_results, data.frame(
      Fish.ID = fish, 
      mlnd = mlnd, 
      CVmlnd = CVmlnd, 
      low10 = low10, 
      low10pc = low10pc, 
      quant_10 = quant[1],
      quant_15 = quant[2],
      quant_20 = quant[3],
      quant_25 = quant[4],
      quant_30 = quant[5]
    ))
    
  } else {
    # If there's not enough valid data for this fish, store NA results
    SMR_results <- rbind(SMR_results, data.frame(
      Fish.ID = fish, 
      mlnd = NA, 
      CVmlnd = NA, 
      low10 = NA, 
      low10pc = NA, 
      quant_10 = NA,
      quant_15 = NA,
      quant_20 = NA,
      quant_25 = NA,
      quant_30 = NA
    ))
    warning(paste("Not enough valid data for Fish.ID:", fish))
  }
}

#2) Deciding which to use 
# --------------------
# see visualizations code. 

mean(SMR_results$low10)
mean(SMR_results$quant_10)
mean(SMR_results$low10pc)


#grab the one's i'm thinking of using
SMR_results_subset <- SMR_results[,c("Fish.ID", "low10", "quant_10", "low10pc")]


#combine SMR data with metabolism data 
data_SMR_Lab_metabolism_A <- merge(data_SMR_Lab_metabolism, SMR_results_subset, by = "Fish.ID") 

# boxplot. SMR by treatment
ggplot(data_SMR_Lab_metabolism_A,
       aes(x = Treatment.x,
           y = low10,
           fill = factor(TempMeta))) +
  geom_boxplot(color = "black",
               position = position_dodge(width = 0.8),
               width = 0.65) +
  scale_fill_manual(
    values = c("10" = "#1f78b4",    # blue for 10 °C
               "15" = "#e31a1c"),   # red  for 15 °C
    name   = "Temperature (°C)"
  ) +
  labs(
    x = "Treatment",
    y = "Resting Metabolism mgO2/kg")+
  theme_minimal(base_size = 18) +
  theme(
    # remove all grid lines
    panel.grid = element_blank(),
    # add crisp black x and y axes
    axis.line = element_line(color = "black", size = 1),
    
   # plot.title  = element_text(face = "bold", size = 22, hjust = 0.5),
   # axis.title  = element_text(face = "bold", size = 18),
   # axis.text.x = element_text(angle = 45, hjust = 1, size = 16),
   # axis.text.y = element_text(size = 16),
   # legend.title= element_text(face = "bold", size = 16),
   # legend.text = element_text(size = 14)
  )


########################################################################################
#------------------------------------------
  #----------------------------#
  ###Calculating EPOC - 
  #----------------------------#
  ## updated from Christians SDA function
  ## smooths curve using quantile smoothing spline. 
  
#STEP 1: creating EPOC function
calcEPOC <- function(X, Y, my.smr, tau = 0.1, lambda = 0.5, 
                       O2.time.unit = c("hour", "min")) {
    require(quantreg)
    require(Hmisc)
    # Prepare data
    my.data <- data.frame(time = X, MO2 = Y)
    my.data <- na.omit(my.data)
    
    # Ensure time is numeric
    my.data$time <- as.numeric(my.data$time)
    
    # Fit model using quantile smoothing spline (rqss)
    fitsEPOC <- rqss(MO2 ~ qss(time, lambda = lambda), tau = tau, data = my.data)
    
    # Predict MO2 based on the fitted model
    pred.x <- data.frame(time = seq(min(my.data$time), max(my.data$time), length.out = 1000)) 
    pred <- predict(fitsEPOC, pred.x)
    
    # Create the EPOC data frame
    epoc <- data.frame(time = pred.x$time, pred = pred)
    epoc$net <- epoc$pred - my.smr  # net MO2
    
    # Set the threshold when net MO2 falls within 10% of the SMR value
    threshold <- my.smr * 0.1  # 10% of SMR
    
    # Find the time when net MO2 falls within 10% of SMR
    end <- subset(epoc, abs(net) <= threshold)[1, 1]
    
    # Find end MO2 value
    end.mo2 <- subset(epoc, time == end)$pred
    
    # Label the stages: "epoc" while metabolic rate is above SMR, "post-epoc" after
    duringEPOC <- epoc$net > threshold & epoc$time < end
    epoc$status <- "post-epoc"
    epoc$status[duringEPOC] <- "epoc"
    epoc$status[epoc$time == min(epoc$time)] <- "start"
    
    # Subset the data to just "epoc" and "start" stages
    epoc.short <- subset(epoc, status %in% c("start", "epoc"))
    {
      # Duration of EPOC in hours
      epoc.dur.h <- max(epoc.short$time)  
      
      # Peak MO2 during EPOC
      epoc.peak <- epoc.short[epoc.short$pred == max(epoc.short$pred), 1:3]  
      
      # Compute the area under the curve (net MO2)
      area.net <- trap.rule(epoc.short$time, epoc.short$net)
      
      # Ensure O2.time.unit is either "hour" or "min"
      O2.time.unit <- match.arg(O2.time.unit)
      
      # Convert area to minutes if necessary
      if (O2.time.unit == "min") {
        area.net <- area.net * 60
      }
    }
    
    # Return the results as a list
    the.epoc <- data.frame(
      duration = epoc.dur.h, 
      peak = epoc.peak$pred, 
      magnitude = area.net, 
      end.MO2 = as.numeric(end.mo2), 
      SMR = my.smr, 
      tau = tau, 
      lambda = lambda
    )
    
    # Rename the peak column
    names(the.epoc)[which(names(the.epoc) == "peak.pred" )] <- "peak"
    
    # Output list with EPOC curve and EPOC variables
    rq.out <- list(epoc.fit = epoc, epoc.var = the.epoc)
    
    return(rq.out)
  }

# STEP 2: Loop through each Fish.ID to calculate EPOC

# Initialize an empty dataframe to store EPOC results
epoc_results <- data.frame()

# Initialize an empty dataframe to store predicted recovery curves
predicted_data <- data.frame()

# Loop through each unique Fish.ID
for (fish_id in unique(data_SMR_Lab_metabolism_A$Fish.ID)) {
  
  # Subset data for the current fish
  EPOC_fish_data <- subset(data_SMR_Lab_metabolism_A, Fish.ID == fish_id)
  
  # Extract MO2, SMR, and time for the current fish
  mo2_data <- EPOC_fish_data$O2.Kg.L.h
  smr <- unique(EPOC_fish_data$low10)  # Ensure SMR is unique or take the first value
  
  # Check if SMR is consistently a single value
  if (length(smr) != 1) {
    warning(paste("Inconsistent SMR values for Fish ID:", fish_id))
    smr <- smr[1]  # Use the first value if inconsistency is found
  }
  
  time <- EPOC_fish_data$Cumulative.Duration
  
  # Call the calcEPOC function for the current fish
  epoc <- calcEPOC(time, mo2_data, smr, O2.time.unit = "hour")
  
  # Store the EPOC variables in the epoc_results dataframe
  epoc_results <- rbind(epoc_results, data.frame(Fish.ID = fish_id, EPOC = epoc$epoc.var))
  
  # Add predicted recovery curve data to the predicted_data dataframe
  epoc_fit <- epoc$epoc.fit
  epoc_fit$Fish.ID <- fish_id  # Add Fish.ID to the predicted data
  predicted_data <- rbind(predicted_data, epoc_fit)
}

# Now I have two dataframes:
# - epoc_results: Contains the EPOC variables for each fish.
# - predicted_data: Contains the predicted recovery curve data for each fish.


# Rename the 'time' column to 'Cumulative.Duration'
colnames(predicted_data)[colnames(predicted_data) == "time"] <- "Cumulative.Duration"

### Merge EPOC data, MO2 Data

# identify the new columns in epoc_results (excluding Fish.ID)
epoc_new_cols <- setdiff(names(epoc_results), "Fish.ID")

#merge only those new columns with the main dataframe
all_results <- merge(
  data_SMR_Lab_metabolism_A,
  epoc_results[, c("Fish.ID", epoc_new_cols)],
  by = "Fish.ID",
  all.x = TRUE
)

summary_df <- all_results %>%
  group_by(Fish.ID) %>%
  slice(1) %>%              # Keep only the first row per fish
  ungroup() %>%             # Remove grouping
  select(Fish.ID, EPOC.duration, EPOC.peak, EPOC.magnitude, 
         EPOC.end.MO2, EPOC.SMR, Treatment.x, Trial, DateMeta, 
         TempMeta, Chamber, TL, Weight.g, low10) %>%
  mutate(
    delta_peak = EPOC.peak - low10,    # difference
    ratio      = EPOC.peak / low10     # ratio
  )
  


#######################################################################################
#------------------------------------------
  #----------------------------#
  ###Calculating Energy expenditure
  #----------------------------#
glimpse(summary_df)

print(mean(summary_df$EPOC.duration))

# Constants
oxy_J_mg <- 13.56  # J per mg O2 

# Calculate EPOC in J/g
summary_df <- summary_df %>%
  mutate(
    EPOC.J.g    = (EPOC.magnitude * oxy_J_mg)/1000   # oxy_J_mg ≈ 14.1 J per mg O2  
  )


print(min(summary_df$EPOC.J.g, na.rm=TRUE))
print(max(summary_df$EPOC.J.g, na.rm=TRUE))

      
      
#summarize energy use - calculate 10th, 50th, and 90th percentiles
stress_EPOC.J <- summary_df %>%
  dplyr::summarise(
    EPOC.J.g_p10 = quantile(EPOC.J.g, probs = 0.10, na.rm = TRUE),
    EPOC.J.g_p50 = quantile(EPOC.J.g, probs = 0.50, na.rm = TRUE),
    EPOC.J.g_p90 = quantile(EPOC.J.g, probs = 0.90, na.rm = TRUE)
  )
print(stress_EPOC.J)

getwd()
#save to disk for use in bioenergetics model
saveRDS(stress_EPOC.J, "../4_BT_FB4_model/Inputs/stress_EPOC.J.rds")



