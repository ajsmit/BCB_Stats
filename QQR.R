# Load required package
library(quantreg)

# Simulated dataset (replace with your actual SST and wind stress curl data)
set.seed(123)
n <- 100000
wind_stress_curl <- rnorm(n)
SST <- rnorm(n, mean = wind_stress_curl, sd = 1)

# Define the quantiles you want to analyze
quantiles <- seq(0.1, 0.9, by = 0.2)

# Create a storage structure for results
results <- list()

# Loop over quantiles of the predictor (wind stress curl)
for (q in quantiles) {
  # Compute quantile for the predictor (wind stress curl)
  wind_quantile <- quantile(wind_stress_curl, probs = q)

  # Subset data based on the quantile of the predictor
  data_subset <- subset(data.frame(wind_stress_curl, SST), wind_stress_curl <= wind_quantile)

  # Loop over quantiles of the response (SST)
  for (r in quantiles) {
    # Perform quantile regression on the subset of data,
    # i.e. the conditional quantiles of SST (in quantile bands tau)
    # given the wind stress curl of winds in each quantile band, q
    model <- rq(SST ~ wind_stress_curl, tau = r, data = data_subset)

    # Store results
    results[[paste("Quantile", q, "Tau", r)]] <- summary(model)
  }
}

# Access and interpret the stored results
results[[1]]  # Example output for the first quantile combination


# Load required packages
library(tidyverse)

# initialise empty data frame to store results
coef_data <- data.frame(Tau = numeric(), Theta = numeric(), Slope = numeric())

# loop through the `results` list and extract the slope coefficients
for (q in quantiles) {
  for (r in quantiles) {
    # extract the slope coefficient for each combination of tau and theta
    slope_coef <- results[[paste("Quantile", q, "Tau", r)]]$coefficients[2, 1]  # Adjust based on model output structure

    # append to the data frame
    coef_data <- rbind(coef_data, data.frame(Tau = q, Theta = r, Slope = slope_coef))
  }
}

# create wide and long versions of the slope data
coef_wide <- coef_data %>%
  pivot_wider(names_from = Tau, values_from = Slope, names_prefix = "Tau_")

coef_long <- coef_data %>%
  mutate(Tau = factor(Tau), Theta = factor(Theta))

# plot the slopes!
ggplot(coef_long, aes(x = Tau, y = Theta, fill = Slope)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  theme_minimal() +
  labs(title = "Heatmap of Quantile-on-Quantile Regression Slopes",
       x = "Quantiles of Wind Stress Curl (Tau)",
       y = "Quantiles of SST (Theta)",
       fill = "Slope Coefficient") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
