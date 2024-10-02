library(circular)
library(tidyverse)

# uas -- Eastward Near-Surface Wind
# vas -- Northward Near-Surface Wind

bind <- readRDS("misc/example_wind_data.Rds")
bind <- as_tibble(bind) |>
  rename(ws = wind_spd,
         wd = wind_dir)

bind$wd <- atan2(-bind$u, -bind$v) * (180 / pi)

# The raw angle will be measured counter-clockwise from the positive x-axis (east).
# To convert this to the standard meteorological convention (where 0° indicates north,
# 90° east, 180° south, and 270° west), use the following adjustment:
bind$wd <- (270 - bind$wd) %% 360

bind$wd1 <- circular(bind$wd, units = "degrees", template = "geographics")
bind$wd2 <- circular(bind$wd, units = "degrees", modulo = "2pi")

library(openair)

windRose(bind, type = "month")


# With ggplot2

library(colorspace)

# Bin wind speed data for categorisation
bind <- bind %>%
  mutate(wind_spd_binned = cut(ws,
                               breaks = c(0, 2, 4, 6, 8, 10, 12, 14, Inf),
                               labels = c("0-2", "2-4", "4-6", "6-8", "8-10", "10-12", "12-14", "14+")))

# Aggregate the data by wind direction and wind speed category
agg_data <- bind %>%
  group_by(month, wd = round(wd/10)*10, wind_spd_binned) %>%
  summarise(count = n()) %>%
  ungroup()

# Function to create a wind rose plot with facet_wrap
create_facet_wind_rose <- function(data) {
  ggplot(data, aes(x = factor(wd), y = count, fill = wind_spd_binned)) +
    geom_bar(stat = "identity", width = 1, colour = "grey70", linewidth = 0.1) +
    coord_polar(start = 0) +
    scale_x_discrete(drop = FALSE, breaks = seq(0, 360, by = 30)) +
    scale_fill_discrete_sequential(palette = "Plasma", rev = FALSE,
                                   name = "Wind\nSpeed\n(m/s)") +
    facet_wrap(~ month, ncol = 3, scales = "free") +
    labs(x = "Wind Direction (degrees)", y = "Frequency") +
    theme_minimal() +
    theme(legend.position = "right")
}

facet_wind_rose <- create_facet_wind_rose(agg_data)
facet_wind_rose

# Save the combined plot to a file
ggsave("WindRose_FacetWrap_Colour_Corrected.png", facet_wind_rose, width = 12, height = 16)
