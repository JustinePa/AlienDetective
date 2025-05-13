# functions.R
# This script contains the custom functions used by the other scripts


fetch_gbif_data <- function(species,
                            hasCoordinate = TRUE,
                            continent = "europe",
                            basisOfRecord = c("OBSERVATION", "MACHINE_OBSERVATION", "HUMAN_OBSERVATION", "MATERIAL_SAMPLE", "LIVING_SPECIMEN", "OCCURRENCE"),
                            fields = c("decimalLatitude", "decimalLongitude", "year", "month", "country"),
                            limit = 10000,
                            output_dir) {
  
  data_list <- rgbif::occ_search(scientificName = species,
                                 hasCoordinate = hasCoordinate,
                                 continent = continent,
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


# Main function: calculates both sea route and geodesic distances from every downloaded GBIF occurrence to the species occurrence in question
calculate.distances <- function(gbif_occurrences, latitude, longitude, raster_map, cost_matrix){
  
  if (is.null(gbif_occurrences)) return(list(sea_distances = NULL, geodesic_distances = NULL, error_messages = "Input table is NULL"))
  if (nrow(gbif_occurrences) < 1) return(list(sea_distances = NULL, geodesic_distances = NULL, error_messages = "Input table is has no entries"))
  
  tryCatch({
    # Specify the PROJ4 string for WGS84
    proj4_crs <- sp::CRS("+init=EPSG:4326")
    
    # Create SpatialPoints objects from the coordinates
    query_point <- sp::SpatialPoints(cbind(longitude, latitude), proj4string = proj4_crs)
    ref_points <- sp::SpatialPoints(cbind(gbif_occurrences$longitude, gbif_occurrences$latitude), proj4string = proj4_crs)
    
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
      sea_distances <- sea_distances / 1000
      geodesic_distances <- geodesic_distances / 1000
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

# Make histogram of sea distances
plot.dist.sea <- function(species, location, distances, output_dir) {
  
  # make histograms of distances per species, with filtering on distance limit 40000
  plot <- ggplot(distances, aes(x = x)) +
    geom_histogram(binwidth = 50, fill = "blue", color = "black", boundary = 0) +
    labs(title = "Histogram of Distances", x = "Distance in km", y = "Frequency of species") +
    theme_bw() +
    #scale_fill_brewer(palette = "Set1") +
    ggtitle(paste0("Distribution of ", species, " from ", location)) +
    theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold"), # set title font size, placement with hjust
          plot.margin = margin(0.3, 0.3, 0.4, 0.4, "cm"),
          axis.text = element_text(size = 16),           # Set font size for axis numbers
          axis.title = element_text(size = 20)) +         # Set font size for axis titles
    scale_x_continuous(breaks = seq(0, 3000, by = 250), expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    coord_cartesian(xlim = c(0, 3000)) # Use coord_cartesian for setting limits
  
  if(!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  ggsave(filename = file.path(output_dir, paste0(gsub(" ", "_", species), "_from_", location, ".png")), 
         plot = plot, width = 2000, height = 1200, units = "px", dpi = 300)
  return(plot)
}


# Make combined histogram of sea distances and fly distances
plot.dist.both <- function(species, location, distances, output_dir) {
  
  plot <- ggplot(distances, aes(x = x, fill = type)) +
    geom_histogram(binwidth = 50, color="#e9ecef", alpha=0.6, position = 'identity') +
    theme_bw() +
    #scale_fill_brewer(palette = "Set1") +
    labs(x = "Distance in km", y = "Frequency") +
    ggtitle(paste0("Distribution of ", species, " from ", location)) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 20, face = "bold"), # set title font size, placement
      plot.margin = margin(0.3, 0.3, 0.4, 0.4, "cm"),
      axis.text = element_text(size = 16),  # Set font size for axis numbers
      axis.title = element_text(size = 20), # Set font size for title
      legend.title = element_text(size = 18, face="bold"), # Settings for legend title
      legend.text = element_text(size = 16)) +  # settings for legend text
    scale_x_continuous(breaks = seq(0, 3000, by = 250), expand = c(0, 0)) +  # settings for x axis
    scale_y_continuous(expand = c(0, 0)) +
    # used expand to make sure the axes are on the lines of the axes and not above them floating
    coord_cartesian(xlim = c(0, 3000)) + # Use coord_cartesian for setting limits
    # set legend title and labels
    scale_fill_discrete(
      name = "Distance type",
      breaks = c("sea", "geodesic"),
      labels = c("Sea distance", "Geodesic distance")) +
    labs(x = "Distance in km", y = "Frequency")
  
  if(!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  ggsave(filename = file.path(output_dir, paste0(gsub(" ", "_", species), "_from_", location, "_seadist&geodesic.png")), 
         plot = plot, width = 2400, height = 1200, units = "px", dpi = 300)
  return(plot)
}


# Make histograms of locations
plot.dist.by.country <- function(species, location, distances, output_dir) {
  
  plot <- ggplot(distances, aes(x = x, fill = country)) +
    geom_histogram(binwidth = 50, boundary = 0, position = "stack") +  # adjust the binwidth to personal preference
    labs(title = paste0("Frequencies of Sea distances/country for ", species," in ", location),
         x = "Sea distance in km", y = "Frequency of species") +
    theme_bw() +
    #scale_fill_brewer(palette = "Set1") +  # You can choose a different palette if you like
    theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold"), # set title font size, placement
          plot.margin = margin(0.3, 0.3, 0.4, 0.4, "cm"),
          axis.text = element_text(size = 16),           # Set font size for axis numbers
          axis.title = element_text(size = 20),
          legend.title = element_text(size = 14),   # Increase legend title size
          legend.text = element_text(size = 12),    # Increase legend text size
          legend.key.size = unit(1.5, "lines")) +   # Increase legend key size
    scale_x_continuous(breaks = seq(0, 3000, by = 250), expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    coord_cartesian(xlim = c(0, 3000)) # Use coord_cartesian for setting limits
  
  if(!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  ggsave(filename = file.path(output_dir, paste0(gsub(" ", "_", species), "_from_", location, "_by_country.png")), 
         plot = plot, width = 2400, height = 1200, units = "px", dpi = 300)
  return(plot)
}


# Make year categories (Used by plot.dist.by.year function)
assign_year_category <- function(year, year_categories) {
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
  
  plot <- ggplot(distances, aes(x = x, fill = year_category)) +
    geom_histogram(binwidth = 50, boundary = 0, position = "stack") +  # adjust the binwidth to personal preference
    labs(title = paste0("Frequencies of Sea distances/year for ", species," in ", location),
         x = "Sea distance in km", y = "Frequency of species") +
    theme_bw() +
    scale_fill_brewer(palette = "YlOrRd", na.value = "black") + # You can choose a different palette if you like
    theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold"), # set title font size, placement
          plot.margin = margin(0.3, 0.3, 0.4, 0.4, "cm"),
          axis.text = element_text(size = 16),           # Set font size for axis numbers
          axis.title = element_text(size = 20),
          legend.title = element_text(size = 14),   # set legend title size
          legend.text = element_text(size = 12),    # set legend text size
          legend.key.size = unit(1.5, "lines")) +   # set legend key size
    scale_x_continuous(breaks = seq(0, 3000, by = 250), expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    coord_cartesian(xlim = c(0, 3000)) # Use coord_cartesian for setting limits
  
  if(!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  ggsave(filename = file.path(output_dir, paste0(gsub(" ", "_", species), "_from_", location, "_by_year.png")), 
         plot = plot, width = 2400, height = 1200, units = "px", dpi = 300)
  return(plot)
}
