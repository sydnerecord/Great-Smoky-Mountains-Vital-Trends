#########################
##### Vegetation ----
#########################

# Generic function to subset trees data by cycle
subset_by_cycle <- function(data, start_year, end_year, exclude_years = NULL) {
  subset_data <- data[data$Year >= start_year & data$Year <= end_year, ]
  if (!is.null(exclude_years)) {
    subset_data <- subset_data[!subset_data$Year %in% exclude_years, ]
  }
  return(subset_data)
}