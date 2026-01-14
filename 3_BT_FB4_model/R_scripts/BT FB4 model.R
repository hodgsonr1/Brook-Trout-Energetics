#Brook Trout (BT) bioenergetics model
#Purpose: 
#----  this code is built off a LMB bioenergetics model from JWB. 
#----- Code is updated and changed to fit parameters to BT.
#Authors:Ryan Hodgson (RH) & Jake Brownscombe (JWB)
#last updated: Nov 29 2025
library(tidyverse)
library(patchwork)
library(flextable)
library(officer)


#basic model ----
oxycal <- 13560 #standard in FB4
EDP <- 4000 #prey energy density (consider altering)
param <- read.csv("Inputs/Parameters_official.csv")
bt <- param %>% filter(Species=="Brook Trout (juvenile & adult)") 
bt$ED <- 7516.5 #range 4162 to 10871. So set to middle manually
bt$ED <- as.numeric(bt$ED)

##############################------------
# Check and assess Brook trout FB4 model components ------------
##############################------------

#data 
fb <- merge(data.frame(temp=seq(1,35,1)), data.frame(weight=seq(100,1000,5))) #create grid of possible temperatures and weights 
head(fb)

#metabolism 
#EQ2
Z <- log(bt$RQ) * (bt$RTM - bt$RTO)
Y <- log(bt$RQ) * (bt$RTM - bt$RTO + 2)
X <- (Z^2 * (1 + sqrt(1 + 40/Y))^2) / 400

V <- (bt$RTM - fb$temp) / (bt$RTM - bt$RTO)
fb$ft <- ifelse(fb$temp < bt$RTM, V^X * exp(X * (1 - V)), 1e-6)

Rmax <- bt$RA * fb$weight ^ bt$RB
fb$Met.J <- Rmax * fb$ft * fb$weight * oxycal

#metabolism values
ggplot(fb, aes(temp, ft))+geom_point()+
ggplot(fb, aes(temp, Rmax, col=weight))+geom_point()+
ggplot(fb, aes(temp, Met.J, col=weight))+geom_point()

#consumption
#EQ3 
CG1 <- (1/(bt$CTO-bt$CQ))*log((0.98*(1-bt$CK1))/(bt$CK1*0.02))
L1 <- exp(CG1*(fb$temp-bt$CQ))
KA <- (bt$CK1*L1) / (1 + bt$CK1*(L1-1))
CG2 <- (1/(bt$CTL-bt$CTM))*log((0.98*(1-bt$CK4))/(bt$CK4*0.02))
L2 <- exp(CG2*(bt$CTL-fb$temp))
KB <- (bt$CK4*L2) / (1 + bt$CK4*(L2-1))
fb$ft <- KA * KB #temp dependence curve
fb$Cmax <- bt$CA * (fb$weight) ^bt$CB #max consumption
fb$Cons.p <- 1 #max consumption at 1
fb$C <- fb$Cmax * fb$Cons.p * fb$ft #actual daily consumption.
fb$Cons.g <- fb$C*fb$weight #g of prey eaten
fb$Cons.J <- fb$Cons.g*EDP #joules

#consumption values:
ggplot(fb, aes(temp, ft))+geom_point()+
  ggplot(fb, aes(temp, Cmax, col=weight))+geom_point()+
  ggplot(fb, aes(temp, C, col=weight))+geom_point()+
  ggplot(fb, aes(temp, Cons.J, col=weight))+geom_point()

#wastes
#EG equation 2
# Egestion (fecal losses)
fb$Eg <- bt$FA * fb$temp^(bt$FB) * exp(bt$FG * fb$Cons.p) * fb$Cons.J
# Excretion (urinary losses)
fb$Ex <- bt$UA * fb$temp^(bt$UB) * exp(bt$UG * fb$Cons.p) * (fb$Cons.J - fb$Eg)

#SDA
fb$SDA <- bt$SDA *(fb$Cons.J-fb$Eg) 

#growth
fb$growth.J <- fb$Cons.J-fb$Met.J-fb$SDA-fb$Eg-fb$Ex
fb$growth.g <- fb$growth.J/bt$ED

#outputs
head(fb)
ggplot(fb, aes(temp, Cons.J, col=weight))+geom_point()+scale_color_viridis_c()+theme_bw()+
ggplot(fb, aes(temp, Met.J, col=weight))+geom_point()+scale_color_viridis_c()+theme_bw()+
ggplot(fb, aes(temp, SDA, col=weight))+geom_point()+scale_color_viridis_c()+theme_bw()+
ggplot(fb, aes(temp, Eg, col=weight))+geom_point()+scale_color_viridis_c()+theme_bw()+
ggplot(fb, aes(temp, Ex, col=weight))+geom_point()+scale_color_viridis_c()+theme_bw()+
ggplot(fb, aes(temp, growth.g, col=weight))+geom_point()+scale_color_viridis_c()+theme_bw()
##############################------------
# Running grow function ------------
##############################------------
#grow a fish over period of time and add simulated C&R 

#load data from resp and feeding projects 
resp_df <- readRDS("Inputs/stress_EPOC.J.rds")
cons_df <- readRDS("Inputs/p_24h_tbl.rds")

##############################------------
# Set up Temp Df #
##############################------------
#Clean Temp Data (weights will be set per scenario)
temp_df <- read.csv("Inputs/temp_daily_scott_stringer_welcome.csv", header=TRUE)
#add a day counter for each fish. 
temp_df <- temp_df %>% 
  mutate(date_est = as.Date(date_est)) %>% #correct date formating
  group_by(fish_name) %>% 
  arrange(date_est, .bygroup = TRUE) %>% 
  mutate(day = as.integer(date_est - min(date_est))+1) %>% #set day count 
  ungroup 

#grab average of median temp across all fish.
sim_template <- temp_df %>% 
  group_by(day) %>% 
  summarize(temp = mean(mean(median_temp), na.rm = TRUE),
            .groups = "drop"
  )

##############################------------
# Calculate weight-specific p0 values ----
##############################------------
# Load weight groups with VBGF-predicted growth
weight_groups <- readRDS("Inputs/weight_groups_vbgf.rds")
# Contains: Weight_Class, initial_g, final_g

# Calculate p0 for each weight class using fit.p
weight_groups <- weight_groups %>%
  rowwise() %>%
  mutate(
    p0 = fit.p(
      initial_g = initial_g,
      final_g = final_g,
      sim_template = sim_template,
      meta = bt,
      oxycal = oxycal,
      EDP = EDP
    )
  ) %>%
  ungroup()

# Create p0 lookup table
p0_lookup <- weight_groups %>%
  select(Weight_Class, p0)

cat("\n=== Weight-Specific p0 Values ===\n")
print(p0_lookup)

##############################------------
#combine into 18 stress scenarios (3 EPOC x 6 feeding situations)
# Pivot EPOC data to long format
resp_long <- resp_df %>%
  pivot_longer(cols = everything(),
               names_to = "EPOC_level",
               values_to = "EPOC.J.g") %>%
  mutate(EPOC_level = case_when(
    EPOC_level == "EPOC.J.g_p10" ~ "low",
    EPOC_level == "EPOC.J.g_p50" ~ "med",
    EPOC_level == "EPOC.J.g_p90" ~ "high"
  ))

# Cross join to create all 18 combinations
stress_meta <- cons_df %>%
  crossing(resp_long)
  # Note: p_stress will be calculated after adding weight-specific p0

# Create parameter grid for all 72 scenarios (18 stress scenarios x 4 C&R event counts)
# Add weight-specific p0 values and calculate p_stress
param_grid <- stress_meta %>%
  crossing(n_events = c(0, 1, 5, 10)) %>%
  left_join(p0_lookup, by = "Weight_Class") %>%
  mutate(p_stress = F_after * p0)
head(param_grid)

# Set up simulation parameters
ndays <- 120
y_spacing <- 5  # days between C&R events

# Load functions and run all 72 scenarios
source("R_scripts/functions.R")

# Run all simulations and combine results with weight-specific p0
results <- param_grid %>%
  pmap_dfr(function(...) {
    scenario_row <- tibble(...)
    run_scenario(
      scenario_row = scenario_row,
      sim_template = sim_template,
      meta = bt,
      ndays = ndays,
      p0 = scenario_row$p0,  # Use weight-specific p0
      y_spacing = y_spacing,
      oxycal = oxycal,
      EDP = EDP
    )
  })

# View results
head(results)
glimpse(results)
view(results)

# Summary statistics
summary(results)
glimpse(results)

#CALCULATE GROWTH REDUCTION + ENERGETIC COSTS
outcomes <- results %>%
  group_by(Weight_Class) %>%
  mutate(
    # 1) Baselines within each weight class (n_events == 0 rows)
    baseline_growth_g   = first(growth_g[n_events == 0]),
    baseline_pct_growth = first(percent_growth[n_events == 0]),
    baseline_energy     = first(Net_Energy[n_events == 0]),
    
    # 2) Growth reduction in grams 
    growth_reduction_g =  baseline_growth_g - growth_g,
    
    # 3) Percentage of grams lost compared to baseline 
    growth_reduction_pct = abs(growth_reduction_g/baseline_growth_g)*100,
    
    # 4) Net energy cost (how much net energy is lost vs baseline)
    net_energy_cost = Net_Energy - baseline_energy,
  
    #5) Delta Percent Growth loss (baseline-natural)
    pct_loss = baseline_pct_growth - percent_growth,
  ) %>%
  ungroup()

 glimpse(outcomes)
 view(outcomes)
  
#extract baseline growth for each group
baseline_growth <- outcomes %>%
  filter(n_events == 0) %>%
  group_by(Weight_Class) %>%
  slice(1) %>%                  
  ungroup() %>%
  select(
    Weight_Class,
    initial_weight_g,
    final_weight_g,
    growth_g,
    percent_growth,
    Net_Energy
  )
glimpse(outcomes)

#SUMMARY TABLE DATA
single_clean <- outcomes %>%
  filter(n_events == 1) %>% 
  rename(
    `Weight Class`       = Weight_Class,
    Temperature          = Temperature,
    EPOC                 = EPOC_level,
    `Energy Cost J`      = net_energy_cost,
    `Growth_Reduction %` = growth_reduction_pct
  ) %>%
  arrange(`Weight Class`, EPOC, Temperature) %>% 
  select(
    `Weight Class`,
    Temperature,
    EPOC,
    `Energy Cost J`,
    `Growth_Reduction %`
  ) %>%
  mutate(
    `Energy Cost J`      = round(`Energy Cost J`, 2),
    `Growth_Reduction %` = round(`Growth_Reduction %`, 2)
  )
ft <- flextable(single_clean)
ft <- autofit(ft)

doc <- read_docx()
doc <- body_add_flextable(doc, ft)
print(doc, target = "Graphs/single_captures_grouped.docx")

####### Growth reduction graph 
#average across EPOC 
df_mean <- outcomes %>%
  filter(n_events > 0) %>%
  group_by(n_events, Weight_Class, Temperature) %>%
  summarise(
    pct_loss = mean(pct_loss, na.rm = TRUE),
    .groups = "drop"
  )                    # keep only C&R scenarios

plot <- ggplot(
  df,
  aes(x = Weight_Class,
      y = pct_loss,
      fill = Temperature)
) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
  facet_wrap(
    ~ n_events,
    labeller = labeller(n_events = function(x) paste(x, "C&R Events"))
  ) +
  scale_fill_manual(
    values = c("10" = "steelblue", "15" = "coral"),
    name   = "Temperature (°C)"
  ) +
  labs(
    x     = "Weight Class",
    y     = "Δ Percent Growth"
  ) +
  theme_bw(base_size = 16) +
  theme(
    axis.text.x     = element_text(angle = 45, hjust = 1, face= "bold"),
    legend.position = "bottom",
    strip.background = element_rect(fill = "grey90", colour = "grey50"),
    strip.text      = element_text(face = "bold"),
    axis.text = element_text(face = "bold"),
    legend_text = element_text(face= "bold"),
    legend.title = element_text(face="bold"),
    axis.title = element_text(face="bold")
  )

plot
ggsave(
  filename = "Graphs/GrowthReduction_Percent.png",
  plot = plot,
  width = 6.5,      # inches
  height = 4.5,     # inches
  dpi = 300
)

library(dplyr)
glimpse(outcomes)
summarize(outcomes)
mean <- outcomes %>% 
  filter(n_events == 1) %>% 
  summarize(mean_pct_loss = mean(pct_loss),
            max_pct_loss = max(pct_loss),
            min_pct_loss = min(pct_loss),
            mean_J = mean(net_energy_cost))
print(mean)



