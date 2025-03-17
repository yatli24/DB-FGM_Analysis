library(ggplot2)
library(dplyr)

setwd("DB-FGM_Analysis")

# Set parameters
days <- 1:365
set.seed(42)

yearly_params <- data.frame(
  Year = years,
  Amplitude = 10 + rnorm(length(years), 0, 1.5),
  PhaseShift = 70 + rnorm(length(years), 0, 4),
  YearlyOffset = rnorm(length(years), 0, 2)
)

# Generate SST data with enhanced yearly variations
sst_data <- data.frame()

for (i in 1:15) {
  year <- yearly_params$Year[i]
  amp <- yearly_params$Amplitude[i]
  phase <- yearly_params$PhaseShift[i]
  y_offset <- yearly_params$YearlyOffset[i]
  
  seasonal_cycle <- 290 + y_offset + amp * sin(2 * pi * (days - phase) / 365)
  noise <- rnorm(365, mean = 0, sd = 4)
  sst <- seasonal_cycle + noise
  
  sst_data <- bind_rows(sst_data, 
                        data.frame(Day = days, 
                                   SST = sst,
                                   Year = as.factor(year)))
}

# Assign alternating line types
sst_data <- sst_data %>%
  mutate(LineType = ifelse(as.numeric(Year) %% 2 == 0, "solid", "dashed"))

# Plot
sst_figure <- ggplot(sst_data, aes(x = Day, y = SST, group = Year, 
                     color = Year, linetype = LineType)) +
  geom_line(alpha = 0.7, linewidth = 0.4) +
  labs(x = "Days since March 1st", y = "SST (degrees K)",
       title = "Observed SST Data at Longitude=-160, Latitude=30 (Simulated)") +
  theme_minimal() +
  theme(legend.position = "none")

# save figure
ggsave("figures/sst_figure_3.png", sst_figure, width = 8, height = 6, dpi = 300)



