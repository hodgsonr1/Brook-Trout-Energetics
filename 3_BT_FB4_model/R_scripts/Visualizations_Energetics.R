## Visualizations Energetic ###

#Description: visualize outputs from FB4 simualtions of catch and release scenarios. 

#Instructions: run and load BT4 model and functions code to produce growth simulation outputs. 

#Ryan Hodgson 
#Nov 29 2025

glimpse(outcomes)
#Growth by weight class x temperature x EPOC level x repeated captures
ggplot(outcomes, aes(
  x = n_events,
  y = growth_reduction_pct,
  color = Weight_Class
)) +
  # Points
  geom_point(size = 4, alpha = 0.9) +

  # Connecting lines
  geom_line(
    aes(group = interaction(Weight_Class, EPOC_level, Temperature)),
    linewidth = 1.2,
    alpha = 0.6
  ) +

  # Facets
  facet_grid(
    rows = vars(Temperature),
    cols = vars(EPOC_level),
    labeller = labeller(
      Temperature = function(x) paste0(x, " °C"),
      EPOC_level  = function(x) paste("EPOC:", x)
    )
  ) +

  ] +

  # Theme
  theme_bw(base_size = 16) +   # bigger base font
  theme(
    strip.text      = element_text(size = 18, face = "bold"),
    axis.title      = element_text(size = 20, face = "bold"),
    axis.text       = element_text(size = 16),
    legend.title    = element_text(size = 18, face = "bold"),
    legend.text     = element_text(size = 16),

    panel.border    = element_rect(color = "black", linewidth = 1),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(linewidth = 0.4, linetype = "dashed"),

    plot.title      = element_text(size = 22, face = "bold", hjust = 0.5)
  ) +

  # Labels
  labs(
    x = "Number of C&R Events",
    y = "Relative Growth Reduction (%)",
    color = "Weight Class"
  )

#Faceted Barplots for Growth Reduction %
outcomes_filtered <- outcomes %>% filter(n_events != 0)

ggplot(outcomes_filtered, aes(
  x = factor(n_events),
  y = growth_reduction_pct,
  fill = (Weight_Class)
)) +
  geom_col(position = "dodge", alpha = 0.9) +

  facet_grid(
    rows = vars(Temperature),
    cols = vars(EPOC_level),
    labeller = labeller(
      Temperature = function(x) paste0(x, " °C"),
      EPOC_level  = function(x) paste("EPOC:", x)
    )
  ) +

  theme_bw(base_size = 16) +
  theme(
    strip.text      = element_text(size = 18, face = "bold"),
    axis.title      = element_text(size = 20, face = "bold"),
    axis.text       = element_text(size = 16),
    axis.text.x     = element_text(angle = 45, hjust = 1),
    legend.title    = element_text(size = 18, face = "bold"),
    legend.text     = element_text(size = 16),

    panel.border    = element_rect(color = "black", linewidth = 1),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(linewidth = 0.4, linetype = "dashed")
  ) +

  labs(
    x = "C&R Events",
    y = "Growth Reduction (%)",
    fill = "Weight Class"
  )

#Summary Table: Net Energy Costs






