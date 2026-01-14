##################################################################
# Title: Respirometry Analysis Kenauk 2024


# Purpose:This script transforms and combines respriometry trial data with meta data, and calculates slopes and background rates.
#         Data consists of 24hr respirometry trials that were run on groups of fish subjected to C&R events of either 
#         (HIGH STRESS) - 30s chase + 10s air exposure, (LOW STRESS) 30s chase only, Or (Control), across 3 temperature groups.


# Date created: October 18 2024
# Date last modified: Nov 18 2024 17:02:29 2024 ------------------------------


#
#Author: Ryan Hodgson
#Contributions from: Paul Bznonek, Christian Bihun and Jake Brownscombe
##################################################################






#Load packages and set global environment
#----------------------------#
#List the packages needed for script
list.of.packages <- c('tidyverse', 'plotly', 'rio', 'readxl')
#Identify packages in the last that are not already on the computer
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
#Install packages in "new.packages"
if(length(new.packages)) install.packages(new.packages); rm(list.of.packages); rm(new.packages)


#Load packages
library(tidyverse) #ggplot, dplyr, etc.
library(plotly) #interactive plots
library(rio) #Import lists of files
library(readxl) #Real excel files
library(respR)




#Set global environment
theme_set(theme_classic()) #Set global ggplot theme
options(scipen=999) #Remove scientific notation


#check working directory 
getwd()




#####Read in the Data#########################################----
#-------------------------------------------------------------# 
###Get the Meta data----#
Lab_MetaData <- read_xlsx("01 Good Raw Data/MetaDataV2_BrookTrout_KenaukResp.xlsx") 
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






#### OPTIONAL CHECK TO SEE IF SOURCES Align between the two dataframes ########


#OPTIONAL: check to make sure all source files are in each dataframe
sourcesRaw <- (unique(data_SMR_Lab$source)) 
print(sourcesRaw)


sourcesMeta <- (unique(Lab_MetaData$source))


all(sourcesRaw %in% sourcesMeta) #Check To some files are named incorrectly


#check for mismatched names
#mismatched_sources <- sourcesRaw[!sourcesRaw %in% sourcesMeta]
#print(mismatched_sources)


##############


## Add metadata, noting that the days wont match-up for overnight trials
# ---------------------------_ #######


#combine the two dataframes 
Combined_data_SMR_Lab <- left_join(data_SMR_Lab, Lab_MetaData, 
                          by = c("source","channel"), relationship = "many-to-many")


# remove all whitespaces
Combined_data_SMR_Lab <- Combined_data_SMR_Lab %>%
  mutate(datetime = str_trim(datetime),
         Trial_Start = str_trim(Trial_Start),
         Trial_End = str_trim(Trial_End))


# keep only selected columns and rename 
Combined_data_SMR_Lab <-Combined_data_SMR_Lab %>%
  select(Fish.ID, Trial, TempMeta=Temp.y, DateMeta=Date.y, DateActual=Date.x, datetime,         
         Treatment, BackgroundVsTrial, Chamber, Channel,              
         Trial_Start, Trial_End, Phase, TL, Weight.g,   
         TempActual=Temp.x, O2, period, type2, Barometric, Notes, source)


########### Optional Quality Check ###


#identifying the removed rows.... 
removed_rows <- Combined_data_SMR_Lab %>%
  filter(!(datetime >= Trial_Start & datetime <= Trial_End))


# filter to look at 1 file at a time
PreBlank_6_Removed <- removed_rows %>% 
  filter(source == "01 Good Raw Data/Trial 6 Aug 15 2024/2024-08-15-PreBlankTrial6_raw.csv") #check for one source
#
#check for NAs 
missing_datetime <- PreBlank_6_Removed %>%
  filter(is.na(datetime))
##################################################


### filtering the combined data 
#___________________###
# filter the new file to keep only values within start and end times
Combined_data_SMR_LabV1 <-Combined_data_SMR_Lab %>%
  filter(datetime >= Trial_Start,
         datetime <= Trial_End) 
#OPTIONAL: check to see if all files remain
print(unique(Combined_data_SMR_LabV1$source))


gc() #Free unused memory


#########################################################


#ensure correct formatting of date time
Combined_data_SMR_LabV1 <- Combined_data_SMR_LabV1 %>%
  mutate(datetime = as.POSIXct(datetime, format = "%Y-%m-%d %H:%M:%S"))


#Trim start and end times of each loop by 30s
df_timebuffer <- Combined_data_SMR_LabV1 %>% 
  group_by(Fish.ID, period) %>% 
  summarize(Trial.Start = min(datetime),
            Trial.End = max(datetime),
            Trial.Start.buffer = Trial.Start + 30,
            Trial.End.buffer = Trial.End - 30,
Reference.Temp = first(TempActual)) %>%    # Get the reference temperature
  group_by(Fish.ID) %>% 
  mutate(Rest.Duration =difftime(Trial.Start, min(Trial.Start), units="hours"))


Combined_data_SMR_LabV1 <- left_join(Combined_data_SMR_LabV1, df_timebuffer, 
                          by = c("Fish.ID", "period")) %>% 
  filter(datetime>=Trial.Start.buffer & datetime<=Trial.End.buffer) 


#### Add temperature adjustment here

## Make function co calculate full saturation
oxygen_saturation <- function(T) {
  14.621 - 0.41022 * T + 0.007991 * T^2 - 0.000077774 * T^3
}

# Rename for simplicity
df<-Combined_data_SMR_LabV1

df$full_sat_reference <- oxygen_saturation(df$Reference.Temp) # Apply function to reference temperature 
df$full_sat_actual <- oxygen_saturation(df$TempActual) # And to all temperatures
df$sat_diff<-df$full_sat_reference -  df$full_sat_actual            # Calculate the difference in O2 saturation between the reference and each timepoint
df$O2_adjusted <- df$O2 - df$sat_diff  #Adjust O2

# Display the dataframe with relevant columns
head(df %>% 
       select(datetime, O2, O2_adjusted, period, Fish.ID, BackgroundVsTrial) %>% 
       as.data.frame())

# Plot both 
b <- ggplot(df %>% 
              filter(Fish.ID == "2024-08-22-037" & BackgroundVsTrial == "Trial")) +
  geom_point(aes(datetime, O2, color = "O2"), size = 0.2) +
  geom_point(aes(datetime, full_sat_actual, color = "Full saturation"), size = 0.3) +
  scale_color_manual(values = c("O2" = "blue", "Full saturation" = "red"))
  
ggplotly(b) # Interactive plot

# Plot both 
b <- ggplot(df %>% 
              filter(Fish.ID == "2024-08-22-037" & BackgroundVsTrial == "Trial")) +
  geom_point(aes(datetime, full_sat_reference, color = "Reference saturation"), size = 0.2) +
  geom_point(aes(datetime, full_sat_actual, color = "Full saturation"), size = 0.3) +
  scale_color_manual(values = c("Reference saturation" = "blue", "Full saturation" = "red"))

ggplotly(b) # Interactive plot


#first make a numeric value column for linear model fitting for each fishID and period:
df <- df %>%
  group_by(Fish.ID, period) %>%
  mutate(datetime.n = row_number()) %>%
  ungroup()

#use one test fish
test.fish <- df %>% filter(Fish.ID=="2024-08-22-037" & BackgroundVsTrial=="Trial")
head(test.fish %>% as.data.frame())

test.fish <- test.fish %>% filter(period == "14")

# Check out one period to make sure it's fitting properly:
p<-ggplot(test.fish, aes(datetime.n, O2)) +
  geom_point(color = "blue", size = 2, alpha = 0.4) +
  #geom_smooth(method = "lm", color = "black", linetype= "dashed") + 
  geom_point(aes(datetime.n, O2_adjusted), color = "red", size = 1, alpha = 0.8)

p

# Check out one period to make sure it's fitting properly:
p<-ggplot(test.fish, aes(datetime.n, full_sat)) +
  geom_point(color = "blue", size = 2, alpha = 0.4) +
  #geom_smooth(method = "lm", color = "black", linetype= "dashed") + 
  geom_point(aes(datetime.n, full_sat_reference), color = "red", size = 1, alpha = 0.8)

p

p<-ggplot(test.fish, aes(datetime.n, sat_diff)) +
  geom_point() +
  geom_smooth(method = "lm") + 
  geom_point(aes(datetime.n, ), color = "red", size = 0.1)

p

p<-ggplot(test.fish, aes(datetime.n, TempActual)) +
  geom_point() +
  geom_smooth(method = "lm") + 
  geom_point(aes(datetime.n, ), color = "red", size = 0.1)

p

# Convert datetime to seconds - relative to each trial
df <- df %>% group_by(Fish.ID, Trial) %>% 
  mutate(time_seconds = as.numeric(difftime(datetime, min(datetime), units = "secs"))) %>%  ungroup()
                                                         
# Get slopes and other stuff
summary_df <- df %>%
  group_by(Fish.ID, period) %>%
  do({
    lm_fit <- lm(O2_adjusted ~ time_seconds, data = .)
    data.frame(
      slope = coef(lm_fit)["time_seconds"],
      Time = first(.$time_seconds) / 3600,  # Convert to hours for the first time point
      r_squared = summary(lm_fit)$r.squared, 
      mean_temp = mean(.$TempActual, na.rm = TRUE),
      sd_temp = sd(.$TempActual, na.rm = TRUE)
    )
  }) %>%
  ungroup()

unique(df$Fish.ID)



#this loop goes through each fish ID and period, selects the data and fits a linear model to it, pulling out a range of O2 values,
#as well as the intercept, slope, r2 of the model, which are the main things you need to calculate O2 consumption
loopC <- c()
for (i in min(test.fish$SS.ID.n, na.rm = T):max(test.fish$SS.ID.n, na.rm = T)){
  loopA <- filter(test.fish, SS.ID.n==i)
  for (j in unique(loopA$period)){
    tryCatch({
      loopB <- filter(loopA, period==j)
      lm <- lm(O2~datetime.n, data=loopB)
      outputs <- cbind(Fish.ID=loopB$Fish.ID[1],
                       period=loopB$period[1],
                       O2max=max(loopB$O2),
                       O2min=min(loopB$O2),
                       O2mean=mean(loopB$O2),
                       O2sd=sd(loopB$O2),
                       O2intercept=lm$coefficients['(Intercept)'], 
                       O2slope=-(lm$coefficients['datetime.n']),
                       r2=summary(lm)$r.squared)
      loopC <- as.data.frame(rbind(loopC, outputs))
    }, error=function(e){NULL})
  }}; rm(loopA); rm(loopB); rm(lm); rm(outputs); rm(i); rm(j)


head(loopC)




#JWB end.   





#####Identify background slopes###################################----
#-------------------------------------------------------------# 


#make dataframe for background trials
background_trials <- df %>% 
  filter(grepl("BLANK", Fish.ID)) %>%
  select(Fish.ID, TempActual, DateMeta, DateActual, datetime, O2, Trial)


# add column  and convert to O2 sat


background_trials<- background_trials %>% 
  mutate(
    O2_sat = convert_DO(
      x = O2,            # Dissolved oxygen value in mg/L
      from = "mg/L",          # Converting from mg/L
      to = "%Air",               # Converting to % saturation
      S = 0,          # Salinity in ppt
      t = TempActual,             # Temperature in °C
      P = 1.013            # Pressure in bar (default atmospheric pressure)
    ))


#Plot each trials background rates in the 4 chambers. (01-04) PRE and (05-08) POST
#Note: trials 10-12 only have 6 because 4th chamber was broken and not used
#change filter(Trial == ""). to whichever trial you wish to view. 




#subset to plot each trial individually
a <- background_trials %>% 
  filter(grepl("BLANK", Fish.ID)) %>% 
  filter(Trial == 6) %>% 
  select(Fish.ID, TempActual, DateMeta, DateActual, datetime, O2, O2_sat)




# graph of O2_sat  
ggplot(sample_n(a, 1000), aes(x=datetime, y=O2_sat))+
  geom_point()+
  facet_wrap(~Fish.ID, scales="free")






#Graph of O2 (mg/L) vs Temp 
ggplot(sample_n(a, 1000), aes(x=datetime)) +
  geom_point(aes(y=O2), color="blue") + 
  geom_line(aes(y=O2), color="blue") + 
  geom_point(aes(y=TempActual), color="red") + 
  geom_line(aes(y=TempActual), color="red") + 
  facet_wrap(~Fish.ID, scales="free") +
  scale_y_continuous(
    name = "O2",
    sec.axis = sec_axis(~ ., name="Temperature")
  ) +
  theme(axis.title.y.left = element_text(color = "blue"),
        axis.title.y.right = element_text(color = "red"))


#NOTE: very strong relationship between the two. Need to control for temp.




##################################################################


#NOTE: I have not progressed much past this point
###############################################################
#Make new table with Manually specified background slopes 


T9_Pre_background <- Combined_data_SMR_LabV1 %>% 
  filter(grepl("BLANK01|BLANK02|BlANK03|BLANK04", Fish.ID) & 
           Trial == 9 &
           datetime > as.POSIXct("2024-08-20 15:45:00") &
           datetime < as.POSIXct("2024-08-20 16:05:00"))






T9_Post_background <- Combined_data_SMR_LabV1 %>% 
  filter(grepl("BLANK05|BLANK06|BlANK07|BLANK08", Fish.ID) & 
           Trial == 9 &
           datetime > as.POSIXct("2024-08-22 15:45:00") &
           datetime < as.POSIXct("2024-08-22 16:00:00"))%>%
  mutate(time_numeric = as.numeric(difftime(datetime, min(datetime), units = "secs")))








#####Analyse the RESP data#########################################----
#-------------------------------------------------------------# 


#----------------------------#
##Use calc_rate to grab slopes # CODE NOT WORKING
#----------------------------#


#convert datetime to "numeric"
#Combined_data_SMR_LabV1$datetime_numeric <- as.numeric(Combined_data_SMR_LabV1$datetime - min(Combined_data_SMR_LabV1$datetime))


# Reorder the columns so that datetime_numeric and O2 are the first two
#Combined_data_SMR_LabV1 <- Combined_data_SMR_LabV1 %>%
  #select(datetime_numeric, O2, everything())  # Select datetime_numeric and O2 first, then all other columns




  
#run through to grab slopes from each loop for each trial
#Combined_data_SMR_LabV1 <- Combined_data_SMR_LabV1 %>%
 # group_by(Fish.ID, period) %>%
  #dplyr::mutate(O2_rate = calc_rate(Combined_data_SMR_LabV1,
                             #from = min(datetime_numeric), 
                             #to = max(datetime_numeric),
                             #by = "time"))%>%
  #ungroup()


#Error in `dplyr::mutate()`:
#ℹ In argument: `O2_rate = calc_rate(...)`.
#ℹ In group 1: `Fish.ID = "2024-08-09 - BLANK 01"` and `period = 1`.
#Caused by error:
 # ! `O2_rate` must be a vector, not a <calc_rate> object.
#Run `rlang::last_trace()` to see where the error occurred.


 #########################          
# PAULS METHOD ## I dont really understand what this code does..


#Loop through experiments and get summary lm results
loopC <- c()
for (i in min(data_SMR_Lab$SS.ID.n, na.rm = T):max(data_SMR_Lab$SS.ID.n, na.rm = T)){
  loopA <- filter(data_SMR_Lab, SS.ID.n==i)
  for (j in unique(loopA$period)){
    tryCatch({
      loopB <- filter(loopA, period==j)
      lm <- lm(O2~datetime.n, data=loopB)
      outputs <- cbind(SS.ID=loopB$SS.ID[1],
                       period=loopB$period[1],
                       O2max=max(loopB$O2),
                       O2min=min(loopB$O2),
                       O2mean=mean(loopB$O2),
                       O2sd=sd(loopB$O2),
                       O2intercept=lm$coefficients['(Intercept)'], 
                       O2slope=-(lm$coefficients['datetime.n']),
                       r2=summary(lm)$r.squared)
      loopC <- rbind(loopC, outputs)
    }, error=function(e){NULL})
  }}; rm(loopA); rm(loopB); rm(lm); rm(outputs); rm(i); rm(j)


  
