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
module load R/4.4.1-cpeGNU-23.12
R --version

# Set custom R library path
export R_LIBS_USER=/cfs/klemming/home/p/pagnier/Rlibs

# Set custom temporary directory
export TMPDIR=/cfs/klemming/home/p/pagnier/tmp
mkdir -p $TMPDIR

# Prevent R from assuming RStudio
export RSTUDIO=0

#set directory
cd /cfs/klemming/home/p/pagnier/testAD

# Run the R script with arguments
Rscript /cfs/klemming/home/p/pagnier/testAD/1.Calculations.R /cfs/klemming/home/p/pagnier/testAD/Species_Location.csv /cfs/klemming/home/p/pagnier/testAD/Coordinates.csv /cfs/klemming/home/p/pagnier/testAD/Output_calculations
