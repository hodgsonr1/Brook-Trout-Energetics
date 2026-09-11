#Title: Temperature Input Growth Model

#Description:
#Script to load in telemetry temperature data for stocked hatchery Brook trout from Little Bent Lake

#Author: Ryan Hodgson

#Date created: Aug 24 2026
#Last updated: Aug 24 2026

#load packages
library(tidyverse)
library(lubridate)
library(data.table)  # only used for fread() - this file is 222MB, base read.csv is much slower

#load data
temp <- fread("Inputs/LittleBent_Temperature.csv") %>% as_tibble()
glimpse(temp)

#format (base tidyverse mutate - temp is a tibble here, not a data.table, so := does not apply)
temp <- temp %>%
  mutate(
    time_EDT = as.POSIXct(time_EDT, tz = "America/Toronto"),
    start_time = as.POSIXct(start_time, tz = "America/Toronto")
  )

#fish w usable temperature. n = 8 (C3 excluded - its readings run consistently
#colder than the other 7 control fish, likely a sensor issue - see earlier project notes)
usable_t <- c("C1", "C2", "C4", "C6", "C5", "C7", "C8", "C10")

#pull controls and temperature data only
temp_clean <- temp %>%
  filter(Units == "°C",
         treatment == "control",
         fish_id %in% usable_t) %>%
  select(fish_id, time_EDT, data, start_time, end_time)
glimpse(temp_clean)

#visualize - all 8 selected fish's raw temperature traces
p_temp_traces <- ggplot(temp_clean, aes(x = time_EDT, y = data, color = fish_id)) +
  geom_line(alpha = 0.7) +
  facet_wrap(~fish_id, ncol = 1, scales = "free_x") +
  labs(x = "Date", y = "Body temperature (°C)",
       title = "Control fish temperature traces") +
  theme_classic(base_size = 12) +
  theme(legend.position = "none")
print(p_temp_traces)
ggsave("Graphs/revisions Sept 2026/LittleBent_TempTraces_AllFish.png", p_temp_traces, width = 8, height = 12, dpi = 300)

##############################------------
# Build a clean 30-day daily composite for the growth model ----
##############################------------
#day-index each fish from its own first recorded day, matching the same convention
temp_clean <- temp_clean %>%
  mutate(date_est = as.Date(time_EDT)) %>%
  group_by(fish_id) %>%
  arrange(date_est, .by_group = TRUE) %>%
  mutate(day = as.integer(date_est - min(date_est)) + 1) %>%
  ungroup()

#daily mean per fish first, then average across fish per day restricted to the first 30 days
daily_by_fish <- temp_clean %>%
  filter(day <= 30) %>%
  group_by(fish_id, day) %>%
  summarize(daily_mean = mean(data, na.rm = TRUE), .groups = "drop")

sim_template_littlebent <- daily_by_fish %>%
  group_by(day) %>%
  summarize(
    temp = mean(daily_mean, na.rm = TRUE),
    n_fish = n(),
    .groups = "drop"
  ) %>%
  arrange(day)

cat("\n=== Little Bent 30-day composite: coverage check ===\n")
print(sim_template_littlebent, n = 30)
cat("\nMinimum fish contributing on any single day:", min(sim_template_littlebent$n_fish), "\n")

#visualize the composite against the individual fish traces, for a sanity check
p_composite <- ggplot() +
  geom_line(data = daily_by_fish, aes(x = day, y = daily_mean, group = fish_id),
            color = "grey70", linewidth = 0.4) +
  geom_line(data = sim_template_littlebent, aes(x = day, y = temp),
            color = "#B85A28", linewidth = 1.2) +
  labs(x = "Day", y = "Daily mean temperature (°C)",
       title = "Little Bent 30-day composite (orange) vs. individual fish (grey)") +
  theme_classic(base_size = 12)
print(p_composite)
ggsave("Graphs/revisions Sept 2026/LittleBent_30Day_Composite.png", p_composite, width = 8, height = 5, dpi = 300)

#save clean day/temp series for BT FB4 model.R
write.csv(
  sim_template_littlebent %>% select(day, temp),
  "Inputs/temp_daily_littlebent_30day.csv",
  row.names = FALSE
)
cat("\nSaved: Inputs/temp_daily_littlebent_30day.csv\n")
