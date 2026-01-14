##Respirometry Visualizations
#Ryan Hodgson
#Last udpated: Oct 29 2025

library(ggplot2)
library(dplyr)
library(patchwork)
library(ggsignif)
library(cowplot) 

### PLOTTING RAW O2

b <- data_SMR_Lab_metabolism%>% 
  filter(TempMeta == 15 & Trial == 1)
str(data_SMR_Lab_metabolism)

ggplot(b, aes(x = Cumulative.Duration, y = O2.Kg.L.h )) +
  geom_point(aes(color = Trial), size = 1) +  # Use lines and color by Trial
  facet_wrap(~ Fish.ID, scales = "free_x") +  # One plot per Fish.ID
  labs(
    title = "O2 by Fish ID",
    x = "Time",
    y = "O2.Kg.L.h",
    color = "Trial")+
  theme_classic(base_size = 12)


##########################################################################################  
##### Background Data #################################################----
#-------------------------------------------------------------# 

####
#plotting 15C background trials
##------------------------------------
str(data_background)

#subset to plot each trial individually
a <- data_background%>% 
  filter(TempMeta == 15 & Trial ==5)

#faceted plot
ggplot(a, aes(x = datetime, y = O2_adjusted)) +
  geom_point(aes(color = Trial), size = 1) +  # Use lines and color by Trial
  facet_wrap(~ Fish.ID, scales = "free_x") +  # One plot per Fish.ID
  labs(
    title = "O2 by Fish ID",
    x = "Period",
    y = "O2 Adjusted",
    color = "Trial")+
  theme_classic(base_size = 12)


###########----------------
single_chamber <- a %>%
  filter(Fish.ID=="2024-08-09 - BLANK 01")

ggplot(a, aes(x = datetime, y = O2_adjusted)) +
  geom_point(aes(color = Trial), size = 1) +  # Use lines and color by Trial
  facet_wrap(~ Fish.ID, scales = "free_x") +  # One plot per Fish.ID
  labs(
    title = "O2 Adjusted Over Period by Fish ID",
    x = "Period",
    y = "O2 Adjusted",
    color = "Trial"
  ) +
  theme_minimal(base_size = 15) +
  theme(
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
  #Plot each trials background rates in the 4 chambers. (01-04) PRE and (05-08) POST
  #Note: trials 10-12 only have 6 because 4th chamber was broken and not used
  #change filter(Trial == ""). to whichever trial you wish to view. 
  
 
  #group into pre and post trials
  a <- a %>%
    mutate(BlankPeriod = ifelse(Fish.ID %in% unique(Fish.ID)[1:4], "Pre Blanks", "Post Blanks"))
  
  #create a label to add "chamber" before the chamber value.
  custom_labeller <- labeller(
    Chamber = function(x) paste("Chamber", x)
  )
  
  #Graph of O2_adjusted (mg/L) vs Temp 
  ggplot(sample_n(a, 1000), aes(x=datetime)) +
    geom_point(aes(y=O2_adjusted), color="blue") + 
    geom_line(aes(y=O2_adjusted), color="blue") + 
    geom_point(aes(y=TempActual), color="red") + 
    geom_line(aes(y=TempActual), color="red") + 
    facet_wrap(~Fish.ID, scales="free") +
    scale_y_continuous(
      name = "O2",
      sec.axis = sec_axis(~ ., name="Temperature")
    ) +
    theme(axis.title.y.left = element_text(color = "blue"),
          axis.title.y.right = element_text(color = "red"))
  #graph O2 
  ggplot(sample_n(a, 1000), aes(x = datetime, y = DO_percent_actual, color = BlankPeriod)) +
    geom_point() +
    facet_wrap(~Fish.ID+ BlankPeriod + Chamber, scales = "free", labeller = custom_labeller) +
    scale_color_manual(values = c("Pre Blanks" = "cadetblue", "Post Blanks" = "cadetblue1"))+
    ggtitle("Trial 1 Oxygen Saturation by Chamber")+
    theme(plot.title = element_text(hjust = 0.5, size = 14)) # Customize the title style
  theme_minimal()
  
  ### check out one blank trial for graphing
  # Filter for the specific Fish.ID
  
  # OPTIONAL 
  # Filter for the specific Fish.ID and fit a linear model
 data_SMR_Lab %>%
    filter(Fish.ID == "2024-08-15- BLANK 04") %>%
    ggplot(aes(x = datetime, y = O2_adjusted)) +
    geom_line() +  # Line plot
    geom_smooth(method = "lm", se = FALSE, color = "blue") +  # Linear model fit (without confidence interval)
    labs(title = "O2 Adjusted Data for Fish.ID: 2024-08-22 - BLANK 01 with Linear Fit",
         x = "DateTime", 
         y = "O2 Adjusted") +
    theme_minimal()
  
  # Fit a linear model to the data
  lm_fit <- lm(O2_adjusted ~ datetime, data = background_trials %>%
                 filter(Fish.ID == "2024-08-24-042", period == "4"))
  
  # View the summary of the model to get the R^2 value
  summary(lm_fit)$r.squared 
  
  
  ###### OPTIONAL##############
  #FROM BEFORE
  ##graph of O2_sat ##### 
  ggplot(sample_n(a, 1000), aes(x = datetime, y = O2_sat, color = BlankPeriod)) +
    geom_point() +
    facet_wrap(~Fish.ID+ BlankPeriod + Chamber, scales = "free", labeller = custom_labeller) +
    scale_color_manual(values = c("Pre Blanks" = "cadetblue", "Post Blanks" = "cadetblue1"))+
    ggtitle("Trial 1 Oxygen Saturation by Chamber")+
    theme(plot.title = element_text(hjust = 0.5, size = 14)) # Customize the title style
  theme_minimal()
  
  
##########################################################################################  
##### Body mass #################################################----
#-------------------------------------------------------------# 
  
  ggplot(summary_df, aes(x = Treatment.x, y = Weight.g)) +
    geom_boxplot(fill = "skyblue", color = "black") +
    labs(title = "Body Mass by treatment", x = "Treatment", y = "Body mass (g)") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12), 
          axis.text.y = element_text(size = 12),
          axis.title = element_text(size = 14),
          plot.title = element_text(size = 16, hjust = 0.5))
# 

##########################################################################################  
##### RAW MO2 data #################################################
#-------------------------------------------------------------# 

#Explore a test fish 
testfish <- data_SMR_Lab_metabolism %>%
  filter(Fish.ID == "2025-06-10-003")

#plot MO2 over time
ggplot(data = testfish, aes(x = Cumulative.Duration, y = O2.Kg.L.h)) +
  geom_point()+
  geom_hline(yintercept = 65, color = "red", linetype = "dashed")+
  labs(title = "Metabolic Rate over Period",
       x = "Time", 
       y = "Metabolic rate (O2/Kg/L/hr)") +
  theme_minimal()

### plot the data for each treatment 

print(unique(data_SMR_Lab_metabolism$Treatment)) #check the treatment names

#check sample sizes
print(length(unique(data_SMR_Lab_metabolism$Fish.ID[data_SMR_Lab_metabolism$Treatment.y == "control"])))
print(length(unique(data_SMR_Lab_metabolism$Fish.ID[data_SMR_Lab_metabolism$Treatment.y == "chase"])))
print(length(unique(data_SMR_Lab_metabolism$Fish.ID[data_SMR_Lab_metabolism$Treatment.y == "chase + air"])))

#subset control fish 
control <- data_SMR_Lab_metabolism %>% 
  filter(Treatment.x == "control", TempMeta == 15, SMR.valid == TRUE)

#subset chase + air fish 
air <- data_SMR_Lab_metabolism %>% 
  filter(Treatment.x == "chase + air", TempMeta == 15, SMR.valid == TRUE)

#
chase <- data_SMR_Lab_metabolism %>% 
  filter(Treatment.x == "chase", TempMeta == 15, SMR.valid == TRUE)

# plot each treatment (note: change control, chase, or air)
ggplot(chase, aes(x = Cumulative.Duration, y = O2.Kg.L.h)) +
  geom_line() +  # Line plot for metabolic rate over time
  facet_wrap(~ Fish.ID, scales = "free_y") +  # Facet by Fish.ID with independent y-axis
  labs(x = "Period", y = "Metabolic Rate (O2 consumption)", title = "15C chase Treatment") +  # Label axes
  theme_minimal() +  # Minimal theme
  theme(axis.text.x = element_text(angle = 45, hjust = 1))#+  # Rotate x-axis labels if needed
  #geom_hline(aes(yintercept = low10), linetype = "dashed", color = "red") 

######################################################
############################################################
#---------------------
  #SMR AND MMR
#---------------------
############################################################


#generic boxplot for looking at treatment on various metabolic data metric (MMR, SMR, EPOC etc)
ggplot(summary_df, aes(x = Treatment.x, y = low10)) +
  stat_boxplot(geom ='errorbar', width=0.25)+
  geom_boxplot(fill = "skyblue", color = "black") +
  labs(title = "SMR by Treatment", x = "Treatment", y = "SMR (MO2)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12), 
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 14),
        plot.title = element_text(hjust = 0.5, size = 16),
        panel.grid.major.y = element_blank(),  # Removes horizontal grid lines
        panel.grid.minor.y = element_blank(),  # Removes minor horizontal grid lines
        axis.line.y = element_line(color = "grey", size = 0.5),# Adds line on the y-axis
        axis.line.x = element_line(color = "grey", size = 0.5))

#---------------------
############################################################
#EPOC METRICS
#---------------------
############################################################

# Reusable function to create standardized EPOC boxplots
create_epoc_boxplot <- function(data, y_var, y_label, base_size = 24) {
  ggplot(data, aes(x = Treatment.x, y = .data[[y_var]], fill = factor(TempMeta))) +
    geom_boxplot(color = "black", position = position_dodge(width = 0.8),
                 linewidth = 1.0) +  # Thicker boxplot lines
    labs(
      x = "Treatment",
      y = y_label,
      fill = "Temperature (°C)"
    ) +
    scale_fill_manual(
      values = c("10" = "#1f78b4", "15" = "#e31a1c"),
      name = "Temperature (°C)"
    ) +
    theme_classic(base_size = base_size) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(linewidth = 1.5, colour = "black"),  # Thicker axis lines
      axis.text = element_text(face = "bold", colour = "black"),    # Bold axis text
      axis.title = element_text(face = "bold", colour = "black"),   # Bold axis titles
      axis.ticks = element_line(linewidth = 1.0)                    # Thicker tick marks
    )
}

# Create individual plots using the function
smr   <- create_epoc_boxplot(summary_df, "low10", "SMR mgO2/kg/h")
peak  <- create_epoc_boxplot(summary_df, "EPOC.peak", " Peak mgO2/kg/h")
mag   <- create_epoc_boxplot(summary_df, "EPOC.magnitude", "EPOC mgO2/kg")
dur   <- create_epoc_boxplot(summary_df, "EPOC.duration", "EPOC duration hrs")
ratio <- create_epoc_boxplot(summary_df, "ratio", "Factorial scope")
delta <- create_epoc_boxplot(summary_df, "delta_peak",  " Scope mgO2/kg/h")

#-------#-------#-------#-------#-------#-------#-------#-------

# ---add a right-side vertical bracket + optional label ---
add_side_bracket <- function(p, label = NULL,
                             x_offset = 0.55,   # how far to the right of last x level
                             tick_w   = 0.08,   # width of top/bottom ticks (increased)
                             frac_low = 0.02,   # bottom of bracket as fraction of y-span
                             frac_high= 0.08,   # top padding below ymax (fraction of y-span)
                             line_size = 1.2) { # Thicker bracket lines
  gb <- ggplot_build(p)
  pp <- gb$layout$panel_params[[1]]

  # y-range (robust to ggplot2 versions)
  y_rng <- if (!is.null(pp$y.range)) pp$y.range else pp$y$range
  y_span <- diff(range(y_rng))
  y0 <- y_rng[1] + frac_low  * y_span        # bottom of bracket
  y1 <- y_rng[2] - frac_high * y_span        # top of bracket

  # number of x categories to place the bar at the far right
  n_x <- if (!is.null(pp$x$breaks)) length(pp$x$breaks) else {
    if (!is.null(pp$x.range)) length(pp$x.range) else length(pp$x$range)
  }
  x_bar <- n_x + x_offset

  p +
    # stem
    geom_segment(inherit.aes = FALSE,
                 aes(x = x_bar, xend = x_bar, y = y0, yend = y1),
                 linewidth = line_size, colour = "black") +
    # ticks
    geom_segment(inherit.aes = FALSE,
                 aes(x = x_bar - tick_w/2, xend = x_bar + tick_w/2, y = y0, yend = y0),
                 linewidth = line_size, colour = "black") +
    geom_segment(inherit.aes = FALSE,
                 aes(x = x_bar - tick_w/2, xend = x_bar + tick_w/2, y = y1, yend = y1),
                 linewidth = line_size, colour = "black") +
    { if (!is.null(label))
      annotate("text", x = x_bar + 0.15, y = (y0 + y1)/2,
               label = label, angle = 90, size = 5, fontface = "bold")  # Larger, bold text
      else NULL } +
    # avoid clipping & give right margin room
    coord_cartesian(clip = "off") +
    theme(plot.margin = margin(6, 30, 6, 6))  # More right margin for larger text
}

#-------#-------#-------#-------#-------#-------#-------#-------

# --- add horizontal bracket above treatment group for within-treatment comparisons ---
add_horizontal_bracket <- function(p,
                                  treatment_x_pos = 2,      # chase is 2nd treatment
                                  label = NULL,
                                  dodge_width = 0.8,
                                  y_offset_frac = 0.15,     # position from top
                                  tick_h = 0.02,            # tick height
                                  line_size = 0.5) {        # match existing brackets

  # Extract plot parameters
  gb <- ggplot_build(p)
  pp <- gb$layout$panel_params[[1]]

  # Get y-range
  y_rng <- if (!is.null(pp$y.range)) pp$y.range else pp$y$range
  y_span <- diff(range(y_rng))

  # Calculate bracket y-position (near top of plot)
  y_bracket <- y_rng[2] - (y_offset_frac * y_span)
  tick_height <- tick_h * y_span

  # Calculate x-positions for dodged boxes
  # With dodge_width = 0.8 and 2 temperature groups:
  dodge_offset <- dodge_width / 4  # = 0.2
  x_left <- treatment_x_pos - dodge_offset   # 10°C box at 1.8
  x_right <- treatment_x_pos + dodge_offset  # 15°C box at 2.2

  # Build bracket
  p_out <- p +
    # Horizontal line
    geom_segment(inherit.aes = FALSE,
                 aes(x = x_left, xend = x_right,
                     y = y_bracket, yend = y_bracket),
                 linewidth = line_size, colour = "black") +
    # Left tick
    geom_segment(inherit.aes = FALSE,
                 aes(x = x_left, xend = x_left,
                     y = y_bracket, yend = y_bracket + tick_height),
                 linewidth = line_size, colour = "black") +
    # Right tick
    geom_segment(inherit.aes = FALSE,
                 aes(x = x_right, xend = x_right,
                     y = y_bracket, yend = y_bracket + tick_height),
                 linewidth = line_size, colour = "black")

  # Add label if provided
  if (!is.null(label)) {
    y_label <- y_bracket + tick_height + (0.04 * y_span)  # Increased spacing
    x_label <- (x_left + x_right) / 2

    p_out <- p_out +
      annotate("text", x = x_label, y = y_label,
               label = label, size = 5, fontface = "bold")
  }

  # Prevent clipping
  p_out <- p_out + coord_cartesian(clip = "off")

  return(p_out)
}

# --- add brackets to each subplot (use your own p-values/labels) ---
smr_b   <- add_side_bracket(smr,   label = "p < 0.001")
peak_b  <- add_side_bracket(peak,  label = "p < 0.002")
mag_b   <- mag    # No vertical bracket
# Add asterisk and p-value for temperature × chase interaction
dur_b <- dur +
  # Asterisk positioned right above the boxes
  annotate("text", x = 2, y = Inf,
           label = "*",
           size = 8, fontface = "bold", vjust = 2.5) +
  # P-value positioned above the asterisk
  annotate("text", x = 2, y = Inf,
           label = "p = 0.05",
           size = 5, fontface = "bold", vjust = 1) +
  coord_cartesian(clip = "off")
ratio_b <- ratio  # No vertical bracket
delta_b <- delta  # No vertical bracket

# --- your original patchwork, now using the bracketed plots ---
## 1) Build the 2×2×2 grid
grid_6 <- (smr_b | peak_b) / (mag_b | dur_b) / (ratio_b | delta_b)

## 2) Optimized theme for single page viewing with larger, bolder text
base_sz <- 14  # Increased from 9 for better visibility
grid_6 <- grid_6 &
  theme_classic(base_size = base_sz) &
  theme(
    axis.title.x = element_blank(),                                      # we'll add one shared x-title
    axis.title.y = element_text(size = base_sz, face = "bold"),         # Bold, larger y-axis titles
    axis.text.x  = element_text(size = base_sz - 1, face = "bold"),     # Bold x-axis text
    axis.text.y  = element_text(size = base_sz - 1, face = "bold"),     # Bold y-axis text
    axis.line    = element_line(linewidth = 1.2, colour = "black"),     # Thicker axis lines
    axis.ticks   = element_line(linewidth = 0.8),                       # Thicker tick marks
    plot.title   = element_text(size = base_sz + 2, face = "bold"),
    legend.title = element_text(size = base_sz, face = "bold"),         # Bold legend title
    legend.text  = element_text(size = base_sz - 1, face = "bold"),     # Bold legend text
    panel.spacing= unit(1.2, "lines"),                                  # More spacing between panels
    plot.margin  = margin(5, 8, 5, 5)                                   # Optimized margins
  )

## 3) Collect legends and place at bottom (NO tags here)
grid_6 <- grid_6 + plot_layout(guides = "collect") &
  theme(
    legend.position = "bottom",
    legend.key.size = unit(0.8, "lines"),
    legend.margin   = margin(t = 4, b = 4),
    legend.box.margin = margin(t = 6)
  )

## 4) Shared x-axis title for categorical x (Treatment)
x_shared <- cowplot::ggdraw() +
  cowplot::draw_label("Treatment", fontface = "bold", size = base_sz + 4)

## 5) Final patchwork with shared x label row + tag only A–F (7th is blank)
patch <- (grid_6 / x_shared + plot_layout(heights = c(1, 1, 1, 0.08))) +
  plot_annotation(tag_levels = list(c("A","B","C","D","E","F",""))) &
  theme(
    plot.tag = element_text(size = base_sz + 6, face = "bold"),
    plot.tag.position = c(0.02, 0.98)
  )

patch

# --- export for 8.5 × 11 in page ---
ggsave("6_EPOC_Metrics_.png", plot = patch,
       width = 11, height = 8.5, units = "in", dpi = 300)

######################################################
############################################################
#---------------------
  #EPOC DATA VISUALIZATIONS #NO MODEL
 # ---------------------
############################################################

#subset control fish n= 11
control_1 <- all_results %>% 
  filter(Treatment.x == "control", TempMeta == 15, SMR.valid == TRUE)

#subset chase + air fish n=11
air_1 <- all_results %>% 
  filter(Treatment.x == "chase + air", SMR.valid == TRUE)

#subset chase fish n=10
chase_1 <- all_results %>% 
  filter(Treatment.x == "chase", SMR.valid == TRUE)


# Create a faceted ggplot highlight AUC 
ggplot(control_1, aes(x = Cumulative.Duration, y = O2.Kg.L.h)) +
  # Plot the raw MO2 data
  geom_line() +
  # Add a horizontal line for SMR (low10)
  geom_hline(aes(yintercept = low10), linetype = "dashed", color = "red") +
  # Add a vertical line for the endpoint
  geom_vline(aes(xintercept = EPOC.duration), linetype = "dotted", color = "blue") +
  # Highlight the EPOC area under the curve (shaded area between curve and SMR line before endpoint)
  geom_ribbon(data = control_1 %>% 
                filter(Cumulative.Duration <= EPOC.duration), 
              aes(ymin = low10, ymax = O2.Kg.L.h), 
              fill = "lightblue", alpha = 0.5) +
  # Facet by Fish.ID
  facet_wrap(~ Fish.ID, scales = "free_y") +
  # Add labels for weight in the top-left corner of each facet
  geom_text(aes(x = 0, y = Inf, label = paste("Weight:", Weight.g, "g")), 
            vjust = 1.5, hjust = -0.1, size = 3, color = "black") +
  # Add labels
  labs(
    x = "Time (hours)", 
    y = "MO2 (mg O2 / kg * L * hour)", 
    title = "EPOC and Raw MO2 for Each Fish"
  ) +
  # Theme adjustments
  theme_minimal() +
  theme(
    strip.text = element_text(size = 10),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

######################################################
############################################################
#---------------------
#EPOC VISUALIZATIONS WITH MODEL
# ---------------------
############################################################

library(stringr)
############################################################
############################################################
############################################################

# Merge data for SMR and EPOC duration from both dataframes
merged_data <- merge(data_SMR_Lab_metabolism_A, epoc_results, by = "Fish.ID", all.x = TRUE)


# Plot MO2 recovery profiles with modelled recovery curves, SMR, and EPOC.duration
ggplot() +
  # Original MO2 data (black line), still filtered by TempMeta = 15
  geom_line(
    data = merged_data %>% filter(TempMeta == 15),
    aes(x = Cumulative.Duration, y = O2.Kg.L.h, group = Fish.ID),
    color = "black"
  ) +
  # Modelled recovery curve (red line), only Fish.IDs starting with "2025"
  geom_line(
    data = predicted_data %>% filter(str_starts(Fish.ID, "2025")),
    aes(x = Cumulative.Duration, y = pred, group = Fish.ID),
    color = "red"
  ) +
  # Horizontal dotted line for SMR (blue), TempMeta = 15
  geom_hline(
    data = merged_data %>% filter(TempMeta == 15),
    aes(yintercept = low10),
    linetype = "dotted", color = "blue"
  ) +
  # Vertical dashed line for EPOC duration (green), TempMeta = 15
  geom_vline(
    data = merged_data %>% filter(TempMeta == 15),
    aes(xintercept = EPOC.duration),
    linetype = "dashed", color = "green"
  ) +
  facet_wrap(~ Fish.ID, scales = "free_y") +
  labs(
    x = "Time (hours)",
    y = "MO2 (mgO2/Kg/L/h)", 
    title = "MO2 Recovery Profile and Modelled Curve (15C)"
  ) +
  theme_minimal() +
  theme(
    strip.text = element_text(size = 12, face = "bold"),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5)
  )

#######

# Single fish EPOC recovery plot
# Fish ID: 2024-08-04_001

# Filter data for the specific fish
fish_id <- "2024-08-04_001"
single_fish_raw <- merged_data %>% filter(Fish.ID == fish_id)
single_fish_pred <- predicted_data %>% filter(Fish.ID == fish_id)

# Get key metrics for annotations
smr_value <- unique(single_fish_raw$low10)
epoc_dur <- unique(single_fish_raw$EPOC.duration)
fish_weight <- unique(single_fish_raw$Weight.g)
fish_temp <- unique(single_fish_raw$TempMeta)

# Create the plot
p_single <- ggplot() +
  # Raw MO2 data (black line)
  geom_line(data = single_fish_raw,
            aes(x = Cumulative.Duration, y = O2.Kg.L.h),
            color = "black", linewidth = 1.5) +
  # Modelled recovery curve (red line)
  geom_line(data = single_fish_pred,
            aes(x = Cumulative.Duration, y = pred),
            color = "red", linewidth = 1.5, linetype = "solid") +
  # SMR horizontal line (blue dashed)
  geom_hline(yintercept = smr_value,
             linetype = "dashed", color = "blue", linewidth = 1.2) +
  # EPOC duration vertical line (green dashed)
  geom_vline(xintercept = epoc_dur,
             linetype = "dashed", color = "darkgreen", linewidth = 1.2) +
  # Labels and formatting
  labs(
    title = paste("EPOC Recovery Profile"),
    subtitle = paste("Weight:", fish_weight, "g | Temperature:", fish_temp, "°C"),
    x = "Time (hours)",
    y = expression(bold("MO"[2]*" (mgO"[2]*"/kg/h)"))
  ) +
  # Theme for presentation
  theme_classic(base_size = 20) +
  theme(
    plot.title = element_text(size = 24, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 18, hjust = 0.5),
    axis.title = element_text(size = 22, face = "bold"),
    axis.text = element_text(size = 18, face = "bold", color = "black"),
    axis.line = element_line(linewidth = 1.5, color = "black"),
    axis.ticks = element_line(linewidth = 1.2),
    panel.background = element_rect(fill = "white"),
    plot.margin = margin(20, 20, 20, 20)
  )

# Display the plot
print(p_single)

# Optional: Save the plot
ggsave(paste0("EPOC_Recovery_", fish_id, ".png"),
       plot = p_single,
       width = 12, height = 8, units = "in", dpi = 300)


###########################################################################
### Relationship of Weight on MO2 ###
########################################################################
########################################################################
######################################################

################## 


# ### Plotting weight as a function of MMR/SMR
ggplot(summary_df, aes(x = log10(Weight.g), y = log10(MMR))) +  
  geom_point() +  # Scatter plot
  geom_smooth(method = "lm", se = FALSE, color = "blue") +  # Line of best fit (linear model)
  labs(
    title = "Relationship between Fish Weight and Mass Specific MMR",  # Plot title
    x = "log10(Weight (g))",  # X-axis label
    y = "log10(MMR (mgO2/kg/h)"  # Y-axis label
  ) +
  theme_minimal()  # Optional: minimal theme for a clean look

# Plotting density plot of fish body mass
ggplot(summary_df, aes(x = log10(Weight.g)) +
  geom_density(fill = "skyblue", color = "black", alpha = 0.7) +
  labs(title = "Density Plot of Fish Body Mass", x = "Body Mass (g)", y = "Density") +
  theme_minimal()


######################################################
############################################################
#---------------------
# Relationship between Weight and EPOC.duration for each TempXTreamtent group
#---------------------
############################################################
library(ggplot2)
library(emmeans)

# 1. Get predicted slopes (trends) for each treatment × temperature
emm_trends <- as.data.frame(
  emtrends(dur_lm, var = "Weight.g", ~ Treatment.x | TempMeta)
)

# 2. Also get predicted lines (fitted values) for smoother visualization
# We'll sample across the observed weight range for nice fitted lines
pred_df <- as.data.frame(
  emmeans(dur_lm, ~ Weight.g | Treatment.x * TempMeta,
          at = list(Weight.g = seq(min(summary_df$Weight.g),
                                   max(summary_df$Weight.g),
                                   length.out = 100)))
)

# 3. Plot
ggplot() +
  # Fitted regression lines per treatment × temperature
  geom_line(data = pred_df,
            aes(x = Weight.g, y = emmean,
                color = Treatment.x, linetype = Treatment.x),
            linewidth = 1) +
  # Confidence ribbons around each line
  geom_ribbon(data = pred_df,
              aes(x = Weight.g,
                  ymin = lower.CL, ymax = upper.CL,
                  fill = Treatment.x),
              alpha = 0.2, colour = NA) +
  # Observed data points (optional, if not too noisy)
  geom_point(data = summary_df,
             aes(x = Weight.g, y = EPOC.duration,
                 color = Treatment.x),
             alpha = 0.6, size = 2) +
  # Facet by temperature
  facet_wrap(~ TempMeta, labeller = label_both) +
  scale_color_manual(values = c("control" = "#1b9e77",
                                "chase" = "#d95f02",
                                "chase + air" = "#7570b3")) +
  scale_fill_manual(values = c("control" = "#1b9e77",
                               "chase" = "#d95f02",
                               "chase + air" = "#7570b3")) +
  labs(
    x = "Body mass (g)",
    y = "EPOC duration (h)",
    color = "Treatment",
    fill = "Treatment",
    linetype = "Treatment",
    title = "Relationship between Body Mass 
    vs Recovery Duration"
  ) +
  theme_minimal(base_size = 16) +
  theme(
    psanel.grid = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.6),
    legend.position = "right",
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5)
  )


