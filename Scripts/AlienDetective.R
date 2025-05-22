# AlienDetective.R
# Main script


#############
### SETUP ###
#############

# Define number of cores 
num_cores <- 4
if (!is.numeric(num_cores) || num_cores <= 0 || num_cores != floor(num_cores)) {
  stop("Number of cores must be a whole number!")
}

# Set CRAN mirror for downloading packages on server
# options(repos = c(CRAN = "https://cloud.r-project.org"))

#setwd("~/AlienDetective")
source("Scripts/functions.R")

# Reset graphics settings
graphics.off()

# NB! No packages loaded here, only installed if missing. Better to use explicit namespaces instead [e.g. raster::extract() rather than just extract()].
# That way it's easier to maintain the code and see which packages are actually required as development progresses, and you also avoid clashes between
# package namespaces, making sure that the correct function is always used regardless of which other packages the user has installed and loaded.
cat(">>> [INIT] Checking for required packages...\n")
packages <- c("rgbif", "sf", "sp", "gdistance", "geodist", "raster", "fasterize", "ggplot2", "rnaturalearth", "rnaturalearthdata", "dplyr", "spThin", "foreach", "doParallel")
for (package in packages) {
  if(!requireNamespace(package, quietly = TRUE)) {
    install.packages(package)
  }
  suppressPackageStartupMessages(library(package, character.only = TRUE))
}

# If installation fails for a package because it is "not available for your version of R", try installing from source
# install.packages(package, pkgType = "source")
# Install the rnaturalearthhires package for a higher-resolution vector map
# install.packages("devtools")
# remotes::install_github("ropensci/rnaturalearthhires")

cat(">>> [INIT] Reading input data...\n")
args <- commandArgs(trailingOnly = TRUE)
species_location_path <- args[1]  # First argument: path to species file
location_coordinates_path <- args[2]  # Second argument: path to coordinates file
rasterized_path <- args[3]  # Third argument: path to rasterized sf-formatted world map (RDS file)
cost_matrix_path <- # Fourth argument: path to transition matrix (RDS file)
output_dir <- args[5]  # Fifth argument: output directory

# Set file paths
if (length(args) == 0) {
  species_location_path <- file.path("Input", "Species_Location_NIS.csv")
  location_coordinates_path <- file.path("Input", "Coordinates_NIS.csv")
  land_polygons_path <- file.path("Input/land_polygons", "land_polygons.shp")
  rasterized_path <- file.path("Input", "rasterized_land_polygons.rds")
  cost_matrix_path <- file.path("Input", "cost_matrix.rds")
  output_dir <- "Output"
}

# if (!file.exists(land_polygons_path) && !file.exists(rasterized_path)) {
#   stop("Path to either an sf-formatted rasterized rds file or a polygon vector shape file must be provided.")
# }

# Read species-location presence/absence matrix
species_location <- read.csv(species_location_path, sep = ";")
# If there are more than one row per species, remove all but the first row for each species
species_location <- species_location[!duplicated(species_location[1]),]
# Read table of coordinates for every location name (ObservatoryID)
location_coordinates <- read.csv(location_coordinates_path, sep = ";")

# INSERT LIST OF NATIVE SPECIES TO REMOVE NATIVE SPECIES FROM DF LIST

# # Subselect species to run the script for (optional). Can also be used to exclude species, e.g. known natives, by negating the which function
species_subset <- c("Potamopyrgus antipodarum")
species_location <- species_location[which(species_location$Specieslist %in% species_subset),]
#species_location <- species_location[c(2, 10, 57),] # Or subset a few species to try at random

required_columns <- c("decimalLatitude", "decimalLongitude", "year", "month", "country")

#########################
### MAP CONFIGURATION ###
#########################

# Load rasterized world map if it exists, otherwise load custom vector shapefile and rasterize it
cat(">>> [MAP] Loading world map...\n")

if(file.exists(rasterized_path)) {
  r <- readRDS(rasterized_path)
} else {
  # Read vector map as sf object
  #land_polygons <- sf::st_read(land_polygons_path)
  land_polygons <- rnaturalearth::ne_countries(scale = "large", returnclass = "sf")
  cat(">>> [MAP] Rasterizing land polygons...\n")
  # Create raster
  r <- raster::raster(raster::extent(-180, 180, -90, 90), crs = sp::CRS("+init=EPSG:4326"), resolution = 0.1)
  # Rasterize vector map using fasterize
  r <- fasterize::fasterize(land_polygons, r, field = NULL, fun = "max")
  # Set sea cells to value 1 and land cells to NA (Opposite of what fasterize outputs)
  r <- raster::calc(r, function(x) ifelse(is.na(x), 1, NA))
  saveRDS(r, rasterized_path)
  rm(land_polygons)
  cat(">>> [MAP] Rasterization done. Saved raster to \"", file.path(getwd(), rasterized_path), "\"\n")
}

if (file.exists(cost_matrix_path)) {
  cat(">>> [MAP] Loading cost matrix...\n")
  cost_matrix <- readRDS(cost_matrix_path)
} else {
  cat(">>> [MAP] Generating cost matrix...\n")
  # Create a transition object for adjacent cells
  cost_matrix <- gdistance::transition(r, transitionFunction = mean, directions = 16)
  # Set infinite costs to NA to prevent travel through these cells
  cost_matrix <- gdistance::geoCorrection(cost_matrix, type = "c", scl = FALSE)
  # Save transition matrix
  saveRDS(cost_matrix, file = cost_matrix_path)
  cat(">>> [MAP] Saved cost matrix to \"", file.path(getwd(), cost_matrix_path), "\"\n")
}

####################################
### Check input coordinates file ###
####################################
# Check if input coordinates are in sea, if not, move them to sea
cat(">>> [COORD] Checking if input coordinates are in sea ...\n")
for (i in 1:nrow(location_coordinates)) {
  loc_name <- location_coordinates$Observatory.ID[i]
  longitude <- as.numeric(gsub(",", ".", location_coordinates$Longitude[i]))
  latitude <- as.numeric(gsub(",", ".", location_coordinates$Latitude[i]))  
  cat("Checking", loc_name,": latitude", latitude, ", longitude", longitude, "\n")
  
  if (is_on_land(latitude, longitude)) {
    cat(loc_name, "is on land, searching nearest sea coordinates...\n")
    moved <- move_to_sea(latitude, longitude)
    
    if (is.null(moved)) {
      cat("No valid sea coordinates found\n")
      message(loc_name, " is on land, no valid sea coordinates found")
    } else {
      # Update df with coordinates moved point
      location_coordinates$Longitude[i] <- moved$coords[1]
      location_coordinates$Latitude[i] <- moved$coords[2]
      dist <- round((moved$dist/1000), 2)
      cat("Updated", loc_name, "to", location_coordinates$Latitude[i], ", ", location_coordinates$Longitude[i], "; moved", dist, "km.\n")
    }
  } else {
    cat(loc_name, "is already in sea\n")
  }
  cat("\n")
}
cat(">>> [DONE] All coordinates updated to nearest sea point\n")

#############################
### DISTANCES CALCULATION ###
#############################

# For non-parallel execution -> use "for" loop
# For parallel execution -> use "foreach" loop + parallel setup

# Setup parallelisation
cluster <- makeCluster(num_cores)
registerDoParallel(cluster)

# for (species in species_location[,1]) {
foreach(species = species_location[,1],
        .packages = c("dplyr", "raster", "sp", "gdistance", "geodist")) %dopar% {
  species_dir <- file.path(output_dir, gsub(" ", "_", species))
  gbif_occurrences_file <- file.path(species_dir, paste0(gsub(" ", "_", species), ".csv"))
  if (file.exists(gbif_occurrences_file)) {
    cat(">>> [GBIF] Loading GBIF data for", species, "\n")
    gbif_occurrences <- read.csv(gbif_occurrences_file, header = TRUE)
  } else {
    cat(">>> [GBIF] Fetching GBIF data for", species, "\n")
    gbif_occurrences <- fetch_gbif_data(species, fields = required_columns)
    if (!is.null(gbif_occurrences)) {
      if (!dir.exists(species_dir)) {
        dir.create(species_dir, recursive = TRUE)
      }
      write.csv(gbif_occurrences, file = gbif_occurrences_file, row.names = FALSE)
    } else {
      # Skip species that have no GBIF records
      next
    }
  }
  
  cat(">>> [GBIF] Ensuring GBIF occurrence coordinates are at sea\n")
  unique_coords <- unique(gbif_occurrences[c("latitude", "longitude")])
  unique_coords$latitude_moved <- NA
  unique_coords$longitude_moved <- NA
  unique_coords$dist_moved <- NA
  counter_moved <- 0
  counter_failed <- 0
  for (i in 1:nrow(unique_coords)) {
    if (is_on_land(unique_coords$latitude[i], unique_coords$longitude[i])) {
      moved <- move_to_sea(unique_coords$latitude[i], unique_coords$longitude[i])
      if (!is.null(moved)) {
        counter_moved <- counter_moved + 1
        unique_coords$latitude_moved[i] <-moved$coords[2]
        unique_coords$longitude_moved[i] <- moved$coords[1]
        unique_coords$dist_moved[i] <- round((moved$dist/1000), 2)
      }
      else {
        counter_failed <- counter_failed + 1
      }
    }
  }
  cat(counter_moved, "of", nrow(unique_coords), "coordinate pairs were moved to sea.\n")
  if (counter_failed != 0) {
    cat("Moving to sea failed for", counter_failed, "coordinate pairs\n")
    message("Species ", species, ": moving to sea failed for ", counter_failed, " coordinate pairs.")
  }
  
  for (location in colnames(species_location[,-1])) {
    # Skip locations where the species hasn't been detected, determined by a read number cutoff (default 1 read).
    # Ideally, data from multiple marker genes should have been compiled into a single presence/absence table before, so there should only be 1 or 0.
    if (species_location[which(species_location[,1] == species), location] < 1) {
      next
    }
    
    # Get coordinates for the location of observation
    latitude <- as.numeric(gsub(",", ".", location_coordinates[which(location_coordinates$Observatory.ID == location), "Latitude"]))
    longitude <- as.numeric(gsub(",", ".", location_coordinates[which(location_coordinates$Observatory.ID == location), "Longitude"]))
    if (length(latitude) != 1 || length(longitude) != 1) {
      message("Could not retrieve coordinates for location \"", location, "\"")
      next
    }
    
    # Run the distance calculations
    cat(">>> [DIST] Calculating distances to", species, "occurrences from", location, "\n")
    result <- calculate.distances(data = unique_coords,
                                  latitude = latitude,
                                  longitude = longitude,
                                  raster_map = r,
                                  cost_matrix = cost_matrix)
    
    if (!(is.null(result$seaway) & is.null(result$geodesic))) {
      unique_coords[,paste0(location, "_seaway")] <- result$sea_distances
      unique_coords[,paste0(location, "_geodesic")] <- result$geodesic_distances
    } else {
      unique_coords[,paste0(location, "_seaway")] <- NA
      unique_coords[,paste0(location, "_geodesic")] <- NA
    }
    if (!is.null(result$error_messages)) {
      message("Errors in calculation for species: ", species, "\n")
      for (error in result$error_messages) {
        message(error)
      }
    }
  }
  # Joining gbif_occurrences df & unique_coords df together
  gbif_occurrences <- gbif_occurrences %>%
    dplyr::left_join(unique_coords, by = c("latitude", "longitude"))

    # Save to csv file
  write.csv(gbif_occurrences, file = gbif_occurrences_file, row.names = FALSE)
  cat("\n")
}

cat(">>> [DONE] Finished calculating distances for all species.")

stopCluster(cluster)
