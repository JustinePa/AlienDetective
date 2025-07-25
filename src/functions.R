# functions.R
# This script contains the custom functions used by the other scripts


fetch_gbif_data <- function(species,
                            hasCoordinate = TRUE,
                            polygon = "POLYGON((-30 81, 50 81, 50 12, -30 12, -30 81))",
                            basisOfRecord = c("OBSERVATION", "MACHINE_OBSERVATION", "HUMAN_OBSERVATION", "MATERIAL_SAMPLE", "LIVING_SPECIMEN", "OCCURRENCE"),
                            fields = c("decimalLatitude", "decimalLongitude", "year", "month", "country"),
                            limit = 10000,
                            output_dir) {
  
  data_list <- rgbif::occ_search(scientificName = species,
                                 hasCoordinate = hasCoordinate,
                                 geometry = polygon,
                                 basisOfRecord = basisOfRecord,
                                 fields = fields,
                                 limit = limit
                                 )
  
  res <- data.frame()
  # Loop over each data frame in the list, add potentially missing columns, and concatente into a single data frame
  for (i in seq_along(data_list)) {
    
    if (is.null(data_list[[i]]$data)) {
      #message("No ", names(data_list[i]), " type records for ", species)
      next
    } else if (nrow(data_list[[i]]$data) == 0) {
      #message("No ", names(data_list[i]), " type records for ", species)
      next
    }
    missing_columns <- setdiff(fields, colnames(data_list[[i]]$data))
    if (length(missing_columns) > 0) {
      for (col in missing_columns) {
        data_list[[i]]$data[col] <- NA
      }
    }
    data_list[[i]]$data$basisOfRecord <- names(data_list)[i]
    res <- rbind(res, data_list[[i]]$data)
  }
  
  if (nrow(res) > 0) {
    # Reorder columns by the order in the fields argument
    res <- res[,c(fields, "basisOfRecord")]
    # Rename lat/long columns
    colnames(res)[colnames(res) == "decimalLatitude"] <- "latitude"
    colnames(res)[colnames(res) == "decimalLongitude"] <- "longitude"
    # Remove occurrences where longitude or latitude is NA
    res <- res[!is.na(res$latitude) & !is.na(res$longitude),]
    return(res)
  } else {
    error_message <- paste0("No GBIF records found for species \"", species, "\"")
    message(error_message)
    return(NULL)
  }
}

# Function to check if point is on land (TRUE = land, FALSE = sea)
is_on_land <- function(lat, lon) {
  point <- sp::SpatialPoints(cbind(lon, lat), proj4string = sp::CRS(proj4string(r)))
  return(is.na(raster::extract(r, point)))
}

# Function to move point on land to sea
move_to_sea <- function(lat, lon) {
  # Get transition matrix & all connected cells
  trans_matrix <- gdistance::transitionMatrix(cost_matrix)
  connected_cells <- which(rowSums(trans_matrix != 0) > 0)
  connected_coords <- raster::xyFromCell(r, connected_cells)
  
  # Filter to only retain sea cells
  is_sea <- raster::extract(r, connected_coords) == 1
  sea_coords <- connected_coords[is_sea, , drop = FALSE]
  
  if (nrow(sea_coords) == 0) {
    return(NULL)  # failure signal
  } 
  
  point_coords <- c(lon, lat)
  
  for (radius_km in seq(5, 100, 5)) {
    radius_deg <- radius_km / 111 # rough conversion km to degrees
    
    # Borders filter box
    lon_min <- lon - radius_deg
    lon_max <- lon + radius_deg
    lat_min <- lat - radius_deg
    lat_max <- lat + radius_deg
    
    # Define sea cells within radius
    sea_in_radius <- which (sea_coords[,1] >= lon_min & sea_coords[,1] <= lon_max &
                            sea_coords[,2] >= lat_min & sea_coords[,2] <= lat_max)
    
    if (length(sea_in_radius) > 0) {
      sea_coords_radius <- sea_coords[sea_in_radius, , drop = FALSE]
      dists <- geosphere::distVincentySphere(point_coords, sea_coords_radius)
      
      # return nearest seapoint
      nearest_idx <- which.min(dists)
      dist <- dists[nearest_idx]
      new_coords <- sea_coords_radius[nearest_idx, , drop = FALSE]
      return(list(
        coords = as.vector(new_coords),
        dist = dist
      ))
    }
    # Otherwise, continues with next larger radius
  }
  return(NULL) # when no sea point found
}


# Main function: calculates both sea route and geodesic distances from every downloaded GBIF occurrence to the species occurrence in question
calculate.distances <- function(data, latitude, longitude, raster_map, cost_matrix){
  
  if (is.null(data)) return(list(sea_distances = NULL, geodesic_distances = NULL, error_messages = "Input table is NULL"))
  if (nrow(data) < 1) return(list(sea_distances = NULL, geodesic_distances = NULL, error_messages = "Input table is has no entries"))
  
  tryCatch({
    # Specify the PROJ4 string for WGS84
    proj4_crs <- sp::CRS("+init=EPSG:4326")
    
    # Create SpatialPoints objects from the coordinates
    query_point <- sp::SpatialPoints(cbind(longitude, latitude), proj4string = proj4_crs)
    ref_points <- sp::SpatialPoints(cbind(ifelse(is.na(data$longitude_moved), data$longitude, data$longitude_moved),
                                          ifelse(is.na(data$latitude_moved), data$latitude, data$latitude_moved)),
                                    proj4string = proj4_crs)
    
    
    # Get raster cell values of the GBIF occurrence points (1 for sea, Inf for land)
    cell_values <- raster::extract(raster_map, ref_points)
    # Initialize result vectors
    sea_distances <- rep(NA_real_, length(cell_values))
    geodesic_distances <- rep(NA_real_, length(cell_values))
    # Get indexes of the points that are in the sea
    indexes <- which(cell_values == 1L)
    if (length(indexes) > 0) {
      # Subset points that are in the sea
      ref_points_sea <- ref_points[indexes,]
      # Vectorized sea distance calculation to all GBIF occurrences in the sea
      sea_distances[indexes] <- as.numeric(gdistance::costDistance(cost_matrix, query_point, ref_points_sea)[1,])
      # Convert points to simple table format for use with geodist
      query_point_table <- data.frame(lon = sp::coordinates(query_point)[,1],
                                      lat = sp::coordinates(query_point)[,2])
      ref_points_sea_table <- data.frame(lon = sp::coordinates(ref_points_sea)[,1],
                                         lat = sp::coordinates(ref_points_sea)[,2])
      # Vecotrized geodesic distance calculation to all GBIF occurrences in the sea
      geodesic_distances[indexes] <- as.numeric(geodist::geodist(query_point_table, ref_points_sea_table, measure = "geodesic"))
      # Convert distances to kilometres
      sea_distances <- round(sea_distances / 1000, 0)
      geodesic_distances <- round(geodesic_distances / 1000, 0)
    }
    # Return result
    return(list(sea_distances = sea_distances, geodesic_distances = geodesic_distances, error_messages = NULL))
  }, error = function(e) {
    error_messages <- paste0("An error occurred during distance calculation for ", species, " in ", location, ": ", e$message)
    return(list(sea_distances = NULL, geodesic_distances = NULL, error_messages = error_messages))
  })
}


##########################
### PLOTTING FUNCTIONS ###
##########################

country.final <- function(species, distances, output_dir) {
  #hist_info <- hist(long_sea$x, plot = FALSE)
  #max_count <- max(hist_info$counts)
  max_x <- max(long_sea$x)
  plot <- ggplot(long_sea, aes(x = x, fill = location)) +
    geom_histogram(binwidth = 50, boundary = 0, position = "stack") + # adjust the binwidth to personal preference
    labs(title = paste0("Frequencies of Sea distances for ", species," from all ARMS locations" ),
         x = "Sea distance in km", y = "Frequency of species") +
    theme_bw() +
    #scale_fill_brewer(palette = "Set1") +  # You can choose a different palette if you like
    theme(plot.title = element_text(hjust = 0.5, size = 12, face = "bold"), # set title font size, placement
          plot.margin = margin(0.3, 0.3, 0.4, 0.4, "cm"),
          axis.text = element_text(size = 10),           # Set font size for axis numbers
          axis.title = element_text(size = 20),
          legend.title = element_text(size = 14),   # Increase legend title size
          legend.text = element_text(size = 12),    # Increase legend text size
          legend.key.size = unit(1.5, "lines")) +   # Increase legend key size
    scale_x_continuous(breaks = seq(0, max_x*1.1, by = 250), expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    coord_cartesian(xlim = c(0, max_x*1.1)) # Use coord_cartesian for setting limits
  
  if(!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  ggsave(filename = file.path(output_dir, paste0(gsub(" ", "_", species), "_from_all_ARMS_locations", ".png")), 
         plot = plot, width = 2400, height = 1200, units = "px", dpi = 300)
  return(plot)
}

# Make histogram of sea distances
plot.dist.sea <- function(species, location, distances, output_dir) {
  hist_info <- hist(sea_loc_data$x, plot = FALSE)
  max_count <- max(hist_info$counts)
  max_x <- max(sea_loc_data$x)
  # make histograms of distances per species, with filtering on distance limit 40000
  plot <- ggplot(sea_loc_data, aes(x = x, fill = location)) +
    geom_histogram(binwidth = 50, boundary = 0, position = "stack", alpha = 0.7) +  # default is position = "stack"
    labs(title = paste("Distances for", species), x = "Distance (km)", y = "Count") +
    theme_minimal()+
    #scale_fill_brewer(palette = "Set1") +
    ggtitle(paste0("Distribution of ", species, " from", loc)) +
    theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold"), # set title font size, placement with hjust
          plot.margin = margin(0.3, 0.3, 0.4, 0.4, "cm"),
          axis.text = element_text(size = 10),           # Set font size for axis numbers
          axis.title = element_text(size = 20)) +         # Set font size for axis titles
    scale_x_continuous(breaks = seq(0, max_x*1.1, by = 250), expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    coord_cartesian(xlim = c(0, max_x*1.1), ylim = c(0,max_count*1.1))
  
  
  if(!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  ggsave(filename = file.path(output_dir, paste0(gsub(" ", "_", species), "_from_", location, ".png")), 
         plot = plot, width = 2000, height = 1200, units = "px", dpi = 300)
  return(plot)
}

# Make combined histogram of sea distances and fly distances
plot.dist.both <- function(species, location, distances, output_dir) {
  combined_distances <- rbind(sea_loc_data, geo_loc_data)
  hist_info <- hist(combined_distances$x, plot = FALSE)
  max_count <- max(hist_info$counts)
  max_x <- max(sea_loc_data$x)
  plot <- ggplot(combined_distances, aes(x = x, fill = DistanceType)) +
    geom_histogram(binwidth = 50, color="#e9ecef", alpha=0.6, position = 'identity') +
    theme_bw() +
    #scale_fill_brewer(palette = "Set1") +
    labs(x = "Distance in km", y = "Frequency") +
    ggtitle(paste0("Distribution of ", species, " from sea and fly distances")) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 20, face = "bold"), # set title font size, placement
      plot.margin = margin(0.3, 0.3, 0.4, 0.4, "cm"),
      axis.text = element_text(size = 10),  # Set font size for axis numbers
      axis.title = element_text(size = 16), # Set font size for title
      legend.title = element_text(size = 18, face="bold"), # Settings for legend title
      legend.text = element_text(size = 16)) +  # settings for legend text
    scale_x_continuous(breaks = seq(0, max_x*1.1, by = 250), expand = c(0, 0)) +  # settings for x axis
    scale_y_continuous(expand = c(0, 0)) +
    # used expand to make sure the axes are on the lines of the axes and not above them floating
    coord_cartesian(xlim = c(0, max_x*1.1), ylim = c(0,max_count*1.1)) + # Use coord_cartesian for setting limits
    # set legend title and labels
    labs(x = "Distance in km", y = "Frequency", fill = "DistanceType")
  
  if(!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  ggsave(filename = file.path(output_dir, paste0(gsub(" ", "_", species), "_from_", location, "_seadist&geodesic.png")), 
         plot = plot, width = 2400, height = 1200, units = "px", dpi = 300)
  return(plot)
}


# Make histograms of locations
plot.dist.by.country <- function(species, location, distances, output_dir) {
  hist_info <- hist(sea_loc_data$x, plot = FALSE)
  max_count <- max(hist_info$counts)
  max_x <- max(sea_loc_data$x)
  plot <- ggplot(sea_loc_data, aes(x = x, fill = country)) +
    geom_histogram(binwidth = 50, boundary = 0, position = "stack") +  # adjust the binwidth to personal preference
    labs(title = paste0("Frequencies of Sea distances/country for ", species," in ", location),
         x = "Sea distance in km", y = "Frequency of species") +
    theme_bw() +
    #scale_fill_brewer(palette = "Set1") +  # You can choose a different palette if you like
    theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold"), # set title font size, placement
          plot.margin = margin(0.3, 0.3, 0.4, 0.4, "cm"),
          axis.text = element_text(size = 10),           # Set font size for axis numbers
          axis.title = element_text(size = 20),
          legend.title = element_text(size = 14),   # Increase legend title size
          legend.text = element_text(size = 12),    # Increase legend text size
          legend.key.size = unit(1.5, "lines")) +   # Increase legend key size
    scale_x_continuous(breaks = seq(0, max_x*1.1, by = 250), expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    coord_cartesian(xlim = c(0, max_x*1.1), ylim = c(0,max_count*1.1)) # Use coord_cartesian for setting limits
  
  if(!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  ggsave(filename = file.path(output_dir, paste0(gsub(" ", "_", species), "_from_", location, "_by_country.png")), 
         plot = plot, width = 2400, height = 1200, units = "px", dpi = 300)
  return(plot)
}


# Make year categories (Used by plot.dist.by.year function)
assign_year_category <- function(year) {
  if (is.na(year)) {
    return(NA)   # return NA when year is not present
  }
  for (category in year_categories) {
    range <- as.numeric(unlist(strsplit(category, "-"))) # save years as numeric without "-"
    if (year >= range[1] & year < range[2]) {  # if the year falls into this category
      return(category) # return this category
    }
  }
  return(NA) # If year doesn't fall into any category, return NA
}


# Make histograms of year categories
plot.dist.by.year <- function(species, location, distances, output_dir) {
  hist_info <- hist(sea_loc_data$x, plot = FALSE)
  max_count <- max(hist_info$counts)
  plot <- ggplot(sea_loc_data, aes(x = x, fill = year_category)) +
    geom_histogram(binwidth = 50, boundary = 0, position = "stack") +  # adjust the binwidth to personal preference
    labs(title = paste0("Frequencies of Sea distances/year for ", species," in ", location),
         x = "Sea distance in km", y = "Frequency of species") +
    theme_bw() +
    scale_fill_brewer(palette = "YlOrRd", na.value = "black") + # You can choose a different palette if you like
    theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold"), # set title font size, placement
          plot.margin = margin(0.3, 0.3, 0.4, 0.4, "cm"),
          axis.text = element_text(size = 10),           # Set font size for axis numbers
          axis.title = element_text(size = 20),
          legend.title = element_text(size = 14),   # set legend title size
          legend.text = element_text(size = 12),    # set legend text size
          legend.key.size = unit(1.5, "lines")) +   # set legend key size
    scale_x_continuous(breaks = seq(0, 7500, by = 250), expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    coord_cartesian(xlim = c(0, 7500), ylim = c(0,max_count*1.1)) # Use coord_cartesian for setting limits
  
  if(!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  ggsave(filename = file.path(output_dir, paste0(gsub(" ", "_", species), "_from_", location, "_by_year.png")), 
         plot = plot, width = 2400, height = 1200, units = "px", dpi = 300)
  return(plot)
}

# In src/functions.R

# plot_leaflet_map <- function(df, species_name, output_dir = NULL) {
#   if (!requireNamespace("leaflet", quietly = TRUE)) {
#     install.packages("leaflet")
#   }
#   library(leaflet)
#   
#   # Basic validation
#   if (!("latitude" %in% tolower(names(df))) || !("longitude" %in% tolower(names(df)))) {
#     stop("Dataframe must contain 'latitude' and 'longitude' columns")
#   }
#   
#   lat_col <- grep("latitude", names(df), ignore.case = TRUE, value = TRUE)
#   lon_col <- grep("longitude", names(df), ignore.case = TRUE, value = TRUE)
#   
#   map <- leaflet(df) %>%
#     addTiles() %>%
#     addCircleMarkers(
#       lng = ~get(lon_col), lat = ~get(lat_col),
#       popup = ~paste("Year:", df$year),
#       color = "blue", radius = 3, stroke = FALSE, fillOpacity = 0.6
#     ) %>%
#     addLegend("bottomright", colors = "blue", labels = "Observations",
#               title = species_name)
#   
#   # Optional: save HTML map
#   if (!is.null(output_dir)) {
#     html_file <- file.path(output_dir, paste0(gsub(" ", "_", species_name), "_map.html"))
#     htmlwidgets::saveWidget(map, file = html_file, selfcontained = TRUE)
#   }
#   
#   return(map)
# }

