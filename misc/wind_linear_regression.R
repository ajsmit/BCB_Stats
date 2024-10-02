library(tidyverse)
library(ncdf4)
library(tidync)
library(FNN)
library(sp)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(lwgeom) # For st_split function
library(nngeo)


# Read data from file directory
data <- tidync("/Users/ajsmit/Downloads/ccam_8km_197901_202304.nc")

data$transforms$time$time <- as.POSIXct((data$transforms$time$time)*60, origin = "1979-01-01")

wind_load <- function(df, lon1, lon2, lat1, lat2, time_start, time_end) {
  wind_dat <- df %>%
    hyper_filter(lon = between(lon, lon1, lon2),
                 lat = between(lat, lat1, lat2),
                 time = time >= time_start & time <= time_end) %>%
    hyper_tibble(select_var = c("uas", "vas"), force = TRUE, drop = TRUE) %>%
    dplyr::rename(t = time, u = uas, v = vas) %>%
    dplyr::mutate(Hour = format(t, "%H"),
                  t = format(t, "%Y-%m-%d")) %>%
    dplyr::filter(Hour == "14") %>%
    dplyr::mutate(wind_spd = sqrt(u^2 + v^2),
                  wind_dir_trig_to = atan2(u/wind_spd, v/wind_spd) ,
                  wind_dir_trig_to_degrees = wind_dir_trig_to * 180/pi,
                  wind_dir_trig_from_degrees = wind_dir_trig_to_degrees + 180) %>%
    dplyr::select(-wind_dir_trig_to, -wind_dir_trig_to_degrees, -Hour)

  return(wind_dat)
}

load_wind <- wind_load(data, lon1 = 16, lon2 = 20, lat1 = -34.6, lat2 = -33.8, time_start = "1980-01-01", time_end = "2022-12-31")

buffer_func <- function(buffer) {

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

new_wind <- load_wind %>%
  group_by(lon, lat) %>%
  summarise()

points_in_buff <- function(df1, df2, df3) {
  coordinates(df1) <- ~lon + lat
  coordinates(df2) <- ~lon + lat
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

buff_wind <- points_in_buff(new_wind, buff, load_wind)

buff_wind$t <- as.Date(buff_wind$t)

slope_func <- function(df) {

  u_models <- df %>%
    group_by(lon, lat) %>%
    group_modify(~ broom::tidy(lm(u ~ t, data = .x))) %>%
    filter(term == "t") %>%
    select(lon, lat, estimate, p.value) %>%
    mutate(p.value = round(p.value, 3),
           estimate = estimate*365*10) %>%
    rename(u_slope = estimate, u_p_value = p.value)

  v_models <- df %>%
    group_by(lon, lat) %>%
    group_modify(~ broom::tidy(lm(v ~ t, data = .x))) %>%
    filter(term == "t") %>%
    select(lon, lat, estimate, p.value) %>%
    mutate(p.value = round(p.value, 3),
           estimate = estimate*365*10) %>%
    rename(v_slope = estimate, v_p_value = p.value)

  spd_models <- df %>%
    group_by(lon, lat) %>%
    group_modify(~ broom::tidy(lm(wind_spd ~ t, data = .x))) %>%
    filter(term == "t") %>%
    select(lon, lat, estimate, p.value) %>%
    mutate(p.value = round(p.value, 3),
           estimate = estimate*365*10) %>%
    rename(wind_spd_slope = estimate, wind_spd_p_value = p.value)

  dir_models <- df %>%
    group_by(lon, lat) %>%
    group_modify(~ broom::tidy(lm(wind_dir_trig_from_degrees ~ t, data = .x))) %>%
    filter(term == "t") %>%
    select(lon, lat, estimate, p.value) %>%
    mutate(p.value = round(p.value, 3),
           estimate = estimate*365*10) %>%
    rename(wind_dir_slope = estimate, wind_dir_p_value = p.value)

  slope_data <- u_models %>% ### updated here
    left_join(v_models, by = c("lon", "lat")) %>%
    left_join(spd_models, by = c("lon", "lat")) %>%
    left_join(dir_models, by = c("lon", "lat"))

  return(slope_data)

}

slope <- slope_func(buff_wind)

library(viridis)

ggplot(slope, aes(x = lon, y = lat)) +
  geom_tile(aes(fill = wind_spd_slope)) +
  borders("world", regions = c("South Africa", "Namibia")) +
  coord_fixed(xlim = c(11, 22), ylim = c(-35, -17)) +
  scale_fill_viridis(option = "turbo")





