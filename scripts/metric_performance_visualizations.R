# This script generates the figures and tables for the performance comparisons
# of the models using the RData files generated from analyze_simulation_results.R

## access analyze_simulation_results.R if needed 
# (must be run after running model_performance_paper.R or model_performance_LV.R)
# file.edit('Helper_functions/analyze_simulation_results.R')

library(ggplot2)
library(dplyr)
library(gt)

setwd("DB-FGM_Analysis")

# Define exact runtimes
runtime_data <- data.frame(
  Model = c("DBFGM", "BGGM", "PSFGM"),
  Runtime_Seconds = c(3 * 60 * 60, 188, 242)  # 3 hours in seconds, BGGM = 188s, PSFGM = 242s
)

# create plot
runtime_plot <- ggplot(runtime_data, aes(x = Model, y = Runtime_Seconds, fill = Model)) +
  geom_bar(stat = "identity") +
  labs(title = "Model Runtimes for 50 Replicates",
       x = "Model",
       y = "Runtime (seconds)") +
  theme_bw() +
  theme(legend.position = "none")

# Save the plot
ggsave("figures/model_runtimes.png", plot = runtime_plot, width = 6, height = 4, dpi = 300)


# Add model name to each table
dbfgm_table <- dbfgm_table %>% mutate(model = "DBFGM")
psfgm_table <- psfgm_table %>% mutate(model = "PSFGM")
bggm_table <- bggm_table %>% mutate(model = "BGGM")

# Combine all tables
combined_df <- bind_rows(dbfgm_table, psfgm_table, bggm_table)

final_df <- combined_df %>%
  filter(Statistic == "Mean", Scenario == "s1") %>%
  select(model, Value, Metric) %>%
  rename(mean = Value, metric = Metric)

# Create the plot
plot <- ggplot(final_df, aes(x = "", y = mean, fill = model)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~ metric, scales = "free_y") + 
  labs(title = "Mean Values by Model for Each Metric (Original Paper Simulation)",
       x = 'Metric',
       y = "Mean Value",
       fill = "Model") +
  ggplot2::theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

# Save the plot as a PNG file
ggsave("figures/mean_values_plot_original.png", plot = plot, width = 8, height = 6, dpi = 300)

# Create a table from the final_df
table_data <- final_df %>%
  tidyr::pivot_wider(names_from = model, values_from = mean) %>%
  arrange(factor(metric, levels = c("TPR", "MCC", "FPR")))

# Generate the table visualization
table_plot <- table_data %>%
  gt() %>%
  tab_header(
    title = "Mean Values by Model Over 50 Replicates",
    subtitle = "Original Paper Simulation Results"
  ) %>%
  cols_label(
    metric = "Metric",
    DBFGM = "DBFGM",
    PSFGM = "PSFGM",
    BGGM = "BGGM"
  ) %>%
  fmt_number(
    columns = c(DBFGM, PSFGM, BGGM),
    decimals = 3
  ) %>%
  tab_options(
    table.font.size = px(14),
    column_labels.font.weight = "bold"
  )

# Save the table as an image
gtsave(table_plot, "figures/mean_values_table_original.png")

# Add model name to each table
dbfgm_table_1 <- dbfgm_table_1 %>% mutate(model = "DBFGM")
psfgm_table_1 <- psfgm_table_1 %>% mutate(model = "PSFGM")
bggm_table_1 <- bggm_table_1 %>% mutate(model = "BGGM")

# Combine all tables
combined_df <- bind_rows(dbfgm_table_1, psfgm_table_1, bggm_table_1)

final_df <- combined_df %>%
  filter(Statistic == "Mean", Scenario == "s1") %>%
  select(model, Value, Metric) %>%
  rename(mean = Value, metric = Metric)

# Create the plot
plot <- ggplot(final_df, aes(x = "", y = mean, fill = model)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~ metric, scales = "free_y") +  
  labs(title = "Mean Values by Model for Each Metric (Lake Victoria Basin Simulation)",
       x = 'Metric',
       y = "Mean Value",
       fill = "Model") +
  ggplot2::theme_bw() + 
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

# Save the plot as a PNG file
ggsave("figures/mean_values_plot_LV.png", plot = plot, width = 8, height = 6, dpi = 300)


# Create a table from the final_df
table_data <- final_df %>%
  tidyr::pivot_wider(names_from = model, values_from = mean) %>%
  arrange(factor(metric, levels = c("TPR", "MCC", "FPR")))  # Ensure correct row order

# Generate the table visualization
table_plot <- table_data %>%
  gt() %>%
  tab_header(
    title = "Mean Values by Model Over 50 Replicates",
    subtitle = "Lake Victoria Simulation Results"
  ) %>%
  cols_label(
    metric = "Metric",
    DBFGM = "DBFGM",
    PSFGM = "PSFGM",
    BGGM = "BGGM"
  ) %>%
  fmt_number(
    columns = c(DBFGM, PSFGM, BGGM),
    decimals = 3
  ) %>%
  tab_options(
    table.font.size = px(14),
    column_labels.font.weight = "bold"
  )

# Save the table as an image
gtsave(table_plot, "figures/mean_values_table_LV.png")

