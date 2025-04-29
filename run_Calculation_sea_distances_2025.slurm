#!/bin/bash
#SBATCH --job-name=sea_distance
#SBATCH --output=sea_distance_%j.out
#SBATCH --error=sea_distance_%j.err
#SBATCH --mem=64G                # Memory allocation
#SBATCH -A naiss2024-22-1364
#SBATCH -t 12:00:00
#SBATCH -N 1
#SBATCH -c 4
#SBATCH --partition=main

# Load required modules
module load PDC
module load R/4.4.1-cpeGNU-23.12
R --version

# Set custom R library path
export R_LIBS_USER=/cfs/klemming/home/m/mattgiel/Rlibs

# Set custom temporary directory
export TMPDIR=/cfs/klemming/home/m/mattgiel/tmp
mkdir -p $TMPDIR

# Prevent R from assuming RStudio
export RSTUDIO=0

#set directory
cd /cfs/klemming/home/m/mattgiel/Alien_Detective

# Run the R script with arguments
Rscript /cfs/klemming/home/m/mattgiel/Alien_Detective/1_Calculation_sea_distances_2025.R /cfs/klemming/home/m/mattgiel/Alien_Detective/Input/Species_Location_NIS.csv /cfs/klemming/home/m/mattgiel/Alien_Detective/Input/Coordinates_NIS.csv /cfs/klemming/home/m/mattgiel/Alien_Detective/Output_calculations
