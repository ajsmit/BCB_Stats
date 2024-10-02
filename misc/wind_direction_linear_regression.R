library(circular)


bind <- readRDS("Processed_data/wind_data_14h.Rds")

bind$wind_dir <- atan2(-bind$u, -bind$v) * (180 / pi) + 180
bind$wind_dir <- circular(bind$wind_dir, units = "degrees", template = "geographics")



complete_data <- transform(bind, season = as.yearqtr(as.yearmon(paste(year, month, sep = "-")) + 1/12))


bind <- complete_data %>% 
  filter(season != "1980 Q1" & season != "2023 Q1") %>% 
  filter(str_detect(season, 'Q1')) %>% 
  group_by(season, lon, lat) %>% 
  summarise(u = mean(u),
            v = mean(v),
            AWS = mean(AWS),
            wind_spd = mean(wind_spd),
            wind_dir = mean(wind_dir),
            curl = mean(curl))

rm(complete_data)
gc()



buffer_func <- function(buffer){
  
  # Download the coastline data for the entire world
  world_coastline <- ne_download(scale = "large", type = "coastline", category = "physical", returnclass = "sf")
  
  # Define the bounding box for Africa
  africa_bbox <- st_bbox(c(xmin = -20, ymin = -35, xmax = 52, ymax = 38), crs = st_crs(world_coastline))
  
  # Crop the coastline data to the bounding box for Africa
  africa_coastline <- st_crop(world_coastline, africa_bbox)
  
  # The BCLME region:
  bclme_bbox <- st_bbox(c(xmin = 7.91755, ymin = -36.61979, xmax = 19.788742, ymax = -5.811113), crs = st_crs(world_coastline))
  
  # Crop the coastline data to the bounding box
  bclme_coastline <- st_crop(africa_coastline, bclme_bbox)
  
  # Create a buffer of 50 nautical miles (1 nautical mile = 1852 meters)
  buffer_50nm <- st_make_valid(st_union(st_buffer(bclme_coastline, dist = buffer * 1852)))
  
  # Ensure coastline extends beyond buffer: Extend the coastline line if necessary
  # (I added/subtracted 3 degrees to the bbox to get a longer coastline)
  # Create a bbox slightly larger than the BCLME extent to get a longer coastline
  # that will insersect the buffer and thus can be used for splitting it lengthwise
  extended_bbox <- st_bbox(c(xmin = 7.91755 - 3,
                             ymin = -36.61979 - 3,
                             xmax = 19.788742 + 3,
                             ymax = -5.811113 + 3),
                           crs = st_crs(world_coastline))
  
  # Crop the coastline data to the outer bbox; this is the coastline that
  # will be used for splitting the buffer
  extended_coastline <- st_crop(africa_coastline, extended_bbox)
  
  # Convert the outer coastline to a single LINESTRING object
  extended_coastline_line <- st_union(st_cast(extended_coastline, "LINESTRING"))
  extended_coastline_line <- st_make_valid(extended_coastline_line)
  #extended_coastline_line_extended <- st_segmentize(extended_coastline_line, dfMaxLength = 1000000) # Adjust the max length as needed
  
  # Split the buffer using the extended coastline line
  split_buffers <- st_split(buffer_50nm, extended_coastline_line)
  
  # Extract the resulting polygons
  split_buffers_sf <- st_collection_extract(split_buffers)
  
  # Separate into two objects
  offshore_buffer <- split_buffers_sf[1, drop = FALSE]
  inland_buffer <- split_buffers_sf[2, drop = FALSE]
  
  # Ensure buffers are in the original CRS
  offshore_buffer <- st_transform(offshore_buffer, crs = st_crs(world_coastline))
  inland_buffer <- st_transform(inland_buffer, crs = st_crs(world_coastline))
  
  # Create simulated gridded data to test the buffer splitting
  # Define the bounding box for the region
  xmin <- 7.91755
  ymin <- -36.61979
  xmax <- 19.788742
  ymax <- -5.811113
  
  # Generate a grid of points within the bounding box
  lon <- seq(xmin, xmax, length.out = 200)
  lat <- seq(ymin, ymax, length.out = 400)
  grid <- expand.grid(lon = lon, lat = lat)
  
  # Add a simulated temperature column
  set.seed(42)  # For reproducibility
  grid$temperature <- runif(nrow(grid), min = 10, max = 30)
  
  # Convert the data frame to an sf object
  grid_sf <- st_as_sf(grid, coords = c("lon", "lat"), crs = st_crs(world_coastline))
  
  # Extract points within the inland buffer
  points_in_inland_buffer <- st_intersects(grid_sf, inland_buffer, sparse = FALSE)
  inland_data <- grid[apply(points_in_inland_buffer, 1, any), ]
  
  # Extract points within the offshore buffer
  points_in_offshore_buffer <- st_intersects(grid_sf, offshore_buffer, sparse = FALSE)
  offshore_data <- grid[apply(points_in_offshore_buffer, 1, any), ]
  
  
  data <- rbind(inland_data, offshore_data)
  
  return(data)
  
}


buff <- buffer_func(50)


new_wind <- bind %>% 
  group_by(lon, lat) %>% 
  summarise()

points_in_buff <- function(df1, df2, df3){
  
  coordinates(df1) <- ~lon+lat
  
  coordinates(df2) <- ~lon+lat
  
  nn1 = get.knnx(coordinates(df2), coordinates(df1), 1)
  
  
  il = nn1$nn.dist[,1]
  
  df1$dist <- il
  
  new_data <- as.data.frame(df1)
  
  update <- new_data %>% 
    filter(dist <= 0.05)
  
  
  new_update <- update %>% 
    left_join(df3, by = c("lon", "lat")) %>% 
    na.omit() %>% 
    select(-dist)
  
  return(new_update)
  
}

buff_wind <- points_in_buff(new_wind, buff, bind)


grouped_data <- buff_wind %>% group_by(lon, lat)


fit_circular_model <- function(df) {
  # Convert wind direction to radians for trigonometric functions
  df$wind_dir_rad <- conversion.circular(df$wind_dir, units="radians")
  
  # Perform regression separately on sine and cosine components
  sin_model <- lm(sin(wind_dir_rad) ~ season, data = df)
  cos_model <- lm(cos(wind_dir_rad) ~ season, data = df)
  
  # Extract slope and p-value for sine and cosine components
  sin_slope <- coef(summary(sin_model))["season", "Estimate"]
  sin_pvalue <- coef(summary(sin_model))["season", "Pr(>|t|)"]
  
  cos_slope <- coef(summary(cos_model))["season", "Estimate"]
  cos_pvalue <- coef(summary(cos_model))["season", "Pr(>|t|)"]
  
  # Calculate the magnitude of the combined slope and a combined p-value
  combined_slope <- sqrt(sin_slope^2 + cos_slope^2)
  combined_pvalue <- max(sin_pvalue, cos_pvalue)
  
  # Return a dataframe with the results
  return(data.frame(
    longitude = df$lon[1],
    latitude = df$lat[1],
    combined_slope = combined_slope,
    combined_pvalue = combined_pvalue
  ))
}

# Apply the function to each group and combine the results into a single dataframe
final_df <- grouped_data %>%
  do(fit_circular_model(.)) %>%
  ungroup()


  
  library(viridis)
  
  ggplot(final_df, aes(x = lon, y = lat)) +
    geom_tile(aes(fill = combined_slope * 10)) +
    borders("world", regions = c("South Africa", "Namibia")) +
    coord_fixed(xlim = c(11, 21), ylim = c(-35, -17)) +
    scale_fill_viridis(option = "turbo") +
    ggtitle("Wind Direction") +
    xlab("Longitude") +
    ylab("Latitude") +
    theme(legend.position = "bottom") +
    labs(fill = NULL)  

  
 

    