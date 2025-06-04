#!/bin/bash
#SBATCH --job-name=AlienDetective
#SBATCH --output=AD_%j.out
#SBATCH --error=AD_%j.err
#SBATCH --mem=64G                # Memory allocation
#SBATCH -A naiss2024-22-1364
#SBATCH -t 5:00:00
#SBATCH -N 1
#SBATCH -c 4
#SBATCH --partition=main
#SBATCH --mail-type=ALL

# Load required modules
module load PDCOLD/23.12
module load R/4.4.1-cpeGNU-23.12
module load proj
module load gdal
module load udunits
R --version

# Set custom R library path
export R_LIBS_USER=/cfs/klemming/home/m/mattgiel/Rlibs

# Set custom temporary directory
export TMPDIR=/cfs/klemming/home/m/mattgiel/tmp
mkdir -p $TMPDIR

# Prevent R from assuming RStudio
export RSTUDIO=0

#set directory
cd /cfs/klemming/home/m/mattgiel/AlienDetective

# Run the R script with arguments
Rscript /cfs/klemming/home/m/mattgiel/AlienDetective/src/AlienDetective.R
