##### Visualizing Feeding Data ----------------------------
#purpose: 
  #code to visualize and explore feeding latency data from Kenauk 2024 and 2025 latency to feed experiments at 15C and 10C 
    #following various simulated C&R events on brook trout (salvenlinus frontinalis)
  
##Instructions:
  #1) Open feeding exp R project 
  #2) load and open "survival analysis" R script 
  #3) make sure to load in dataframe from above script and clean it
  #4) run this code to visualize latency to feed data

#Author: Ryan Hodgson
#created: Oct 9 2025 
#last modified: Dec 18 2025
###----------------------------

#boxplot 
ggplot(feed_df, aes(x = Treatment,
                    y = Latency_hrs,
                    fill = factor(Temperature..C))) +   # convert to factor for discrete fill
  geom_boxplot(outlier.shape = 21, color = "black") +
  scale_fill_manual(
    values = c("10" = "#1f78b4",   # blue
               "15" = "#e31a1c")   # red
  ) +
  labs(
    x = "Treatment",
    y = "Latency to Feed (hours)",
    fill = "Temperature (°C)"
  ) +
  theme_classic(base_size = 22)


#interaction plot 
summary_df <- feed_df %>%
  group_by(Temperature..C, Treatment) %>%
  summarise(mean_time = mean(Latency_hrs, na.rm = TRUE),
            se = sd(Latency_hrs, na.rm = TRUE)/sqrt(n()), .groups="drop")

ggplot(summary_df, aes(x = Treatment,
                       y = mean_time,
                       group = Temperature..C,
                       color = Temperature..C)) +
  geom_line(position = position_dodge(0.2)) +
  geom_point(size = 3, position = position_dodge(0.2)) +
  scale_fill_manual(
    values = c("10" = "#1f78b4",   # blue
               "15" = "#e31a1c")   # red
  ) +
  geom_errorbar(aes(ymin = mean_time - se,
                    ymax = mean_time + se),
                width = 0.15,
                position = position_dodge(0.2)) +
  labs(y = "Mean Time to Feed (±SE)") +
  theme_classic(base_size = 22)


### Survival Plot Temperature  ####
#instructions: 
#--> load data and cox model from "survival analysis" script.
#--> load the newdata temp_df which keeps weight.kg constant (average), and treatment as control. Since no effect of treatment 
#--> fit survival model with new data
#--> plot 


#construct new df for temperature at average weight
mw <- mean(feed_df$Weight.kg, na.rm = TRUE)
temp_df <- tibble(
  Temperature..C = factor(c("10","15"), levels = levels(feed_df$Temperature..C)),
  Treatment      = factor("control_A", levels = levels(feed_df$Treatment)),
  Weight.kg      = mw
)

#fit model
fit <- survfit(cph, newdata = temp_df)
#plot
pt <- ggsurvplot(
  fit,
  data = feed_df,
  conf.int = TRUE,
  palette = c("#1f78b4", "#e31a1c"),
  legend.title = "Temp (°C)",
  legend.labs  = c("10 °C", "15 °C"),
  xlab = "Time to first feeding (hours)",
  ylab = "Survival",
  ggtheme = theme_classic(base_size = 12) +
    theme(
      # Thick axis lines
      axis.line = element_line(linewidth = 1.5, colour = "black"),
      axis.ticks = element_line(linewidth = 1.2, colour = "black"),
      # Remove all grid lines
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      # Large bold text
      axis.text = element_text(size = 12, face = "bold", colour = "black"),
      axis.title = element_text(size = 12, face = "bold", colour = "black"),
      legend.title = element_text(size = 12, face = "bold"),
      legend.text = element_text(size = 12, face = "bold"),
      # Clean background
      panel.background = element_rect(fill = "white", colour = NA),
      plot.background = element_rect(fill = "white", colour = NA)
    )
)
print(pt)
ggsave("03_Graphs/KM_temp_Surv_Plot.png",width=18, height=18, units="cm", dpi=300)

#Survival Plot Weight 
#instructions: 
#--> load data and cox model from "survival analysis" script.
#--> create the newdata weight_df which grabs 10th and 90th percentile weights, and keeps temp (15C) and treatment as control.
#--> fit survival model with new data
#--> plot 

# Grab 10th and 90th percentiles
weight_q <- quantile(feed_df$Weight.kg, probs = c(0.1, 0.9), na.rm = TRUE)
print(weight_q)
# New weight_q# New df with two weight groups for control fish only
weight_KM <- tibble(
  Temperature..C = factor(c("10"), levels = levels(feed_df$Temperature..C)),
  Treatment      = factor("control_A", levels = levels(feed_df$Treatment)),
  Weight.kg      = weight_q
)

# Fit KM survival curves for each temperature
fit_weight <- survfit(cph, newdata = weight_KM)

pw <- ggsurvplot(
  fit_weight,
  data = feed_df,
  conf.int = TRUE,
  palette = c("#4daf4a", "#984ea3"),
  legend.title = "Weight",
  legend.labs  = c("Small", "Large"),
  xlab = "Time to first feeding (hours)",
  ylab = "Survival",
  ggtheme = theme_classic(base_size = 12) +
    theme(
      # Thick axis lines
      axis.line = element_line(linewidth = 1.5, colour = "black"),
      axis.ticks = element_line(linewidth = 1.2, colour = "black"),
      # Remove all grid lines
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      # Large bold text
      axis.text = element_text(size = 12, face = "bold", colour = "black"),
      axis.title = element_text(size = 12, face = "bold", colour = "black"),
      legend.title = element_text(size = 12, face = "bold"),
      legend.text = element_text(size = 12, face = "bold"),
      # Clean background
      panel.background = element_rect(fill = "white", colour = NA),
      plot.background = element_rect(fill = "white", colour = NA)
    )
)
print(pw)
ggsave("03_Graphs/KM_temp_weight.png",width=18, height=18, units="cm", dpi=300)

#combine plots using patchwork
library(patchwork)
library(cowplot)

p1 <- pw$plot
p2 <- pt$plot

plots <- (p1 | p2) + 
  plot_layout(axis_titles = "collect", guides = "collect") +
  plot_annotation(tag_levels = "A" ) &
  theme(
    plot.tag = element_text(face = "bold"),
    legend.position = "right",
    legend.box = "vertical"
  ) &
  guides(colour = guide_legend(n_col = 1))

print(plots)

ggsave("03_Graphs/SurvPlotPanel.png",width=7.0, height=3.5, units="in", dpi=300)







