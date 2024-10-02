# Load the numDeriv package
library(tidyverse) # For basic data manipulation
library(tidync) # For easily dealing with NetCDF data
library(doMC)
doMC::registerDoMC(cores = 8)

library(numDeriv)

nc_file <- "/Users/ajsmit/Downloads/ccam_8km_197901_202304.nc"

loadnc <- function(file_name, lon1, lon2, lat1, lat2) {
  dat <- tidync(file_name) %>%
    hyper_filter(lon = between(lon, lon1, lon2),
                 lat = between(lat, lat1, lat2)) %>%
    hyper_tibble() %>%
    select(lon, lat, time, vas, uas) %>%
    mutate(time = as.Date(time, origin = "1979-01-01"))
  return(dat)
}

df <- loadnc(nc_file, 12, 19, -34, -32)

# Function to calculate wind stress curl
rho_air <- 1.225  # Air density (kg/m^3)
Cd <- 1.2e-3      # Drag coefficient

# Calculate wind stress components
df <- df |>
  mutate(tau_x = rho_air * Cd * uas * sqrt(uas^2 + vas^2),
         tau_y = rho_air * Cd * vas * sqrt(uas^2 + vas^2))

calculate_wind_stress_curl <- function(df) {
  # Function to calculate the curl at a single time step
  wind_stress_curl <- function(tau_x, tau_y, lon, lat) {
    # Define gradient calculation for tau_y with respect to lon
    d_tau_y_dx <- array(0, dim = dim(tau_y))
    for (i in 1:dim(tau_y)[1]) {
      for (j in 1:dim(tau_y)[2]) {
        d_tau_y_dx[i, j, ] <- grad(function(lon_val) {
          idx <- which.min(abs(lon - lon_val))
          return(tau_y[i, j, idx])
        }, lon)
      }
    }

    # Define gradient calculation for tau_x with respect to lat
    d_tau_x_dy <- array(0, dim = dim(tau_x))
    for (i in 1:dim(tau_x)[1]) {
      for (j in 1:dim(tau_x)[3]) {
        d_tau_x_dy[i, , j] <- grad(function(lat_val) {
          idx <- which.min(abs(lat - lat_val))
          return(tau_x[i, idx, j])
        }, lat)
      }
    }

    # Calculate the curl
    curl <- d_tau_y_dx - d_tau_x_dy
    return(curl)
  }

  # Apply the curl function to each time step
  result <- df %>%
    group_by(time) %>%
    nest() %>%
    mutate(curl = map(data, function(data) {
      lon <- unique(data$lon)
      lat <- unique(data$lat)
      tau_x <- array(data$tau_x, dim = c(length(lon), length(lat)))
      tau_y <- array(data$tau_y, dim = c(length(lon), length(lat)))
      curl <- wind_stress_curl(tau_x, tau_y, lon, lat)
      data$curl <- as.vector(curl)
      return(data)
    })) %>%
    unnest(c(data, curl))

  return(result)
}

result_df <- calculate_wind_stress_curl(df)

df2 <- df |>
  dplyr::filter(time == "1979-12-27")

wind_stress_curl <- function(tau_x, tau_y, lon, lat) {
  # Define gradient calculation for tau_y with respect to lon
  d_tau_y_dx <- array(0, dim = dim(tau_y))
  for (i in 1:dim(tau_y)[1]) {
    for (j in 1:dim(tau_y)[2]) {
      d_tau_y_dx[i, j, ] <- grad(function(lon_val) {
        idx <- which.min(abs(lon - lon_val))
        return(tau_y[i, j, idx])
      }, lon)
    }
  }

  # Define gradient calculation for tau_x with respect to lat
  d_tau_x_dy <- array(0, dim = dim(tau_x))
  for (i in 1:dim(tau_x)[1]) {
    for (j in 1:dim(tau_x)[3]) {
      d_tau_x_dy[i, , j] <- grad(function(lat_val) {
        idx <- which.min(abs(lat - lat_val))
        return(tau_x[i, idx, j])
      }, lat)
    }
  }

  # Calculate the curl
  curl <- d_tau_y_dx - d_tau_x_dy
  return(curl)
}


lon <- unique(df2$lon)
lat <- unique(df2$lat)
tau_x <- array(df2$tau_x, dim = c(length(lon), length(lat)))
tau_y <- array(df2$tau_y, dim = c(length(lon), length(lat)))
curl <- wind_stress_curl(tau_x, tau_y, lon, lat)
df2$curl <- as.vector(curl)
