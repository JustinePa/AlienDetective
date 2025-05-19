# AlienDetective.R
# Main script


#############
### SETUP ###
#############

#setwd("~/AlienDetective")
source("Scripts/functions.R")

# Reset graphics settings
graphics.off()

# NB! No packages loaded here, only installed if missing. Better to use explicit namespaces instead [e.g. raster::extract() rather than just extract()].
# That way it's easier to maintain the code and see which packages are actually required as development progresses, and you also avoid clashes between
# package namespaces, making sure that the correct function is always used regardless of which other packages the user has installed and loaded.
message(">>> [INIT] Checking for required packages...")
packages <- c("rgbif", "sf", "sp", "gdistance", "geodist", "raster", "fasterize", "ggplot2", "rnaturalearth", "rnaturalearthdata", "dplyr", "spThin")
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

message(">>> [INIT] Reading input data...")
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

# Subselect species to run the script for (optional). Can also be used to exclude species, e.g. known natives, by negating the which function
species_subset <- c("Amphibalanus amphitrite")
species_location <- species_location[which(species_location$Specieslist %in% species_subset),]
#species_location <- species_location[c(2, 10, 57),] # Or subset a few species to try at random

required_columns <- c("decimalLatitude", "decimalLongitude", "year", "month", "country")


#########################
### MAP CONFIGURATION ###
#########################

# Load rasterized world map if it exists, otherwise load custom vector shapefile and rasterize it
message(">>> [MAP] Loading world map...")

if(file.exists(rasterized_path)) {
  r <- readRDS(rasterized_path)
} else {
  # Read vector map as sf object
  #land_polygons <- sf::st_read(land_polygons_path)
  land_polygons <- rnaturalearth::ne_countries(scale = "large", returnclass = "sf")
  message(">>> [MAP] Rasterizing land polygons...")
  # Create raster
  r <- raster::raster(raster::extent(-180, 180, -90, 90), crs = sp::CRS("+init=EPSG:4326"), resolution = 0.1)
  # Rasterize vector map using fasterize
  r <- fasterize::fasterize(land_polygons, r, field = NULL, fun = "max")
  # Set sea cells to value 1 and land cells to NA (Opposite of what fasterize outputs)
  r <- raster::calc(r, function(x) ifelse(is.na(x), 1, NA))
  saveRDS(r, rasterized_path)
  rm(land_polygons)
  message(">>> [MAP] Rasterization done. Saved raster to \"", file.path(getwd(), rasterized_path), "\"")
}

if (file.exists(cost_matrix_path)) {
  message(">>> [MAP] Loading cost matrix...")
  cost_matrix <- readRDS(cost_matrix_path)
} else {
  message(">>> [MAP] Generating cost matrix...")
  # Create a transition object for adjacent cells
  cost_matrix <- gdistance::transition(r, transitionFunction = mean, directions = 16)
  # Set infinite costs to NA to prevent travel through these cells
  cost_matrix <- gdistance::geoCorrection(cost_matrix, type = "c", scl = FALSE)
  # Save transition matrix
  saveRDS(cost_matrix, file = cost_matrix_path)
  message(">>> [MAP] Saved cost matrix to \"", file.path(getwd(), cost_matrix_path), "\"")
}

####################################
### Check input coordinates file ###
####################################
# Check if input coordinates are in sea, if not, move them to sea
message(">>> [COORD] Checking if input coordinates are in sea ...")
for (i in 1:nrow(location_coordinates)) {
  loc_name <- location_coordinates$Observatory.ID[i]
  longitude <- as.numeric(gsub(",", ".", location_coordinates$Longitude[i]))
  latitude <- as.numeric(gsub(",", ".", location_coordinates$Latitude[i]))  
  message("Checking ", loc_name, ": latitude ", latitude, ", longitude ", longitude)
  
  # Define the point
  inp_point <- sp::SpatialPoints(cbind(longitude, latitude), proj4string = sp::CRS(proj4string(r)))
  
  if (is_on_land(inp_point)) {
    message(loc_name, " is on land, searching nearest sea coordinates...")
    moved_point <- move_to_sea(inp_point)
    
    if (is.null(moved_point)) {
      message("No valid sea coordinates found for ", loc_name)
    } else {
      # Update df with coordinates moved point
      location_coordinates$Longitude[i] <- sp::coordinates(moved_point)[1]
      location_coordinates$Latitude[i] <- sp::coordinates(moved_point)[2]
      message("Updated ", loc_name, " to latitude ", location_coordinates$Latitude[i], ", longitude ", location_coordinates$Longitude[i])
    }
  } else {
    message(loc_name, " is already in sea")
  }
  message("")
}
message(">>> [DONE] All coordinates updated to nearest sea point")

#############################
### DISTANCES CALCULATION ###
#############################

for (species in species_location[,1]) {
  species_dir <- file.path(output_dir, gsub(" ", "_", species))
  gbif_occurrences_file <- file.path(species_dir, paste0(gsub(" ", "_", species), ".csv"))
  if (file.exists(gbif_occurrences_file)) {
    message(">>> [GBIF] Loading GBIF data for ", species)
    gbif_occurrences <- read.csv(gbif_occurrences_file, header = TRUE)
  } else {
    message(">>> [GBIF] Fetching GBIF data for ", species)
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
  
  message(">>> [GBIF] Ensuring GBIF occurrence coordinates are at sea")
  unique_coords <- unique(gbif_occurrences[c("latitude", "longitude")])
  unique_coords$latitude_moved <- NA
  unique_coords$longitude_moved <- NA
  unique_coords$dist_moved <- NA
  counter <- 0
  for (i in 1:nrow(unique_coords)) {
    if (is_on_land(unique_coords$latitude[i], unique_coords$longitude[i])) {
      moved <- move_to_sea(unique_coords$latitude[i], unique_coords$longitude[i])
      if (!is.null(moved)) {
        counter <- counter + 1
        unique_coords$latitude_moved[i] <-moved$coords[2]
        unique_coords$longitude_moved[i] <- moved$coords[1]
        unique_coords$dist_moved[i] <- round((moved$dist/1000), 2)
      }
    }
  }
  message(counter, " of ", nrow(unique_coords), " coordinate pairs were moved to sea.")

#!!! Need to add code to write to file, but not same file as GBIF data because other dimensions
  #write.csv(gbif_occurrences, file = gbif_occurrences_file, row.names = FALSE)
  
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
      warning("Could not retrieve coordinates for location \"", location, "\"")
      next
    }
    
    # Run the distance calculations
    message(">>> [DIST] Calculating distances to ", species, " occurrences from ", location)
    result <- calculate.distances(gbif_occurrences = thinned_dataframe,
                                  latitude = latitude,
                                  longitude = longitude,
                                  raster_map = r,
                                  cost_matrix = cost_matrix)
    if (!(is.null(result$seaway) & is.null(result$geodesic))) {
      gbif_occurrences[,paste0(location, "_seaway")] <- result$sea_distances
      gbif_occurrences[,paste0(location, "_geodesic")] <- result$geodesic_distances
    } else {
      gbif_occurrences[,paste0(location, "_seaway")] <- NA
      gbif_occurrences[,paste0(location, "_geodesic")] <- NA
    }
    if (!is.null(result$error_messages)) {
      for (error in result$error_messages) {
        message(error)
      }
    }
  }
  
  # Save to csv file
  #write.csv(gbif_occurrences, file = gbif_occurrences_file, row.names = FALSE)
  message("")
}

################
### PLOTTING ###
################

# library("ggplot2")
# required_columns <- c("latitude", "longitude", "year", "month", "country", "basisOfRecord")
# 
# # Iterate over species for which a csv file exists
# for (species in species_location[,1]) {
#   species_ <- gsub(" ", "_", species)
#   species_dir <- file.path(output_dir, species_)
#   if (file.exists(file.path(species_dir, paste0(species_, ".csv")))) {
#     distance_df <- read.csv(file.path(species_dir, paste0(species_, ".csv")))
#   } else {
#     warning("No output directory found for species \"", species, "\". Skipping plotting.")
#     next
#   }
#   
#   # Iterate over locations in that csv file
#   for (location in unique(gsub("_seaway|_geodesic", "", colnames(distance_df[,-which(colnames(distance_df) %in% required_columns)])))) {
#     
#     ### BELOW PART IS UNFINISHED
#     
#     seaway_df <- distance_df[,c(required_columns, paste0(location, "_seaway"))]
#     seaway_df <- reshape(seaway_df,
#                          varying   = vals,           # the columns to stack
#                          v.names   = "value",        # name of the value column
#                          timevar   = "variable",     # name of the column that will hold the former column names
#                          times     = vals,           # the values to put in that new “variable” column
#                          idvar     = "Specieslist",  # keep Specieslist as your ID
#                          direction = "long")
#     
#     both_df <- distance_df[,c(required_columns, paste0(location, "_geodesic"))]
#     both_df <- reshape(both_df,
#                        varying   = vals,           # the columns to stack
#                        v.names   = "value",        # name of the value column
#                        timevar   = "variable",     # name of the column that will hold the former column names
#                        times     = vals,           # the values to put in that new “variable” column
#                        idvar     = "Specieslist",  # keep Specieslist as your ID
#                        direction = "long")
#     
#     # Add distance type column, used by plot.dist.both function
#     seaway_df$type <- "seaway"
#     # Assign year categories
#     year_categories <- c("1965-1985", "1985-1990", "1990-1995",
#                          "1995-2000", "2000-2005", "2005-2010",
#                          "2010-2015", "2015-2020", "2020-2025")
#     seaway_df$year_category <- sapply(seaway_df$year, assign_year_category(year_categories = year_categories))
#     both_df$year_category <- sapply(both_df$year, assign_year_category(year_categories = year_categories))
#     # clean dataframe from rows with Inf in them
#     seaway_df <- seaway_df[is.finite(seaway_df$Koster_seaway), ]
#     # select only distances below 40000km
#     seaway_df <- subset(seaway_df, Koster_seaway < 40000)
#     
#     plot <- plot.dist.sea(
#       species = species,
#       location = location,
#       distances = seaway_df,
#       output_dir = species_dir
#     )
#     
#     plot <- plot.dist.both(
#       species = species,
#       location = location,
#       distances = both_df,
#       output_dir = species_dir
#     )
#     
#     plot <- plot.dist.by.country(
#       species = species,
#       location = location,
#       distances = seaway_df,
#       output_dir = species_dir
#     )
#     
#     plot <- plot.dist.by.year(
#       species = species,
#       location = location,
#       distances = seaway_df,
#       output_dir = species_dir
#     )
#     
#   }
# }
