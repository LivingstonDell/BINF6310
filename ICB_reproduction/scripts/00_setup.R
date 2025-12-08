 # ============================================
# 00_setup.R
# ICB Reproduction Project - Environment Setup
# Author: Trivedi Dhai
# Date: November 2025
# ============================================

# DESCRIPTION:
# This script sets up the R environment with required packages
# Run this first before any other scripts

# ============================================
# LOAD REQUIRED LIBRARIES
# ============================================
library(Seurat)      # Single-cell analysis
library(dplyr)       # Data manipulation
library(ggplot2)     # Plotting

# ============================================
# SET WORKING DIRECTORY
# ============================================
# Uncomment and modify path as needed
# setwd("/scratch/trivedi.dhai/ICB_reproduction")

# ============================================
# DEFINE PATHS
# ============================================
paths <- list(
  raw_data = "/scratch/trivedi.dhai/ICB_reproduction/data/raw/",
  processed_data = "/scratch/trivedi.dhai/ICB_reproduction/data/processed/",
  figures = "/scratch/trivedi.dhai/ICB_reproduction/results/figures/",
  tables = "/scratch/trivedi.dhai/ICB_reproduction/results/tables/",
  scripts = "/scratch/trivedi.dhai/ICB_reproduction/scripts/"
)

# ============================================
# DEFINE COLOR PALETTES
# ============================================
dataset_colors <- c("Melanoma" = "#E41A1C", 
                    "BCC" = "#377EB8", 
                    "Liver" = "#4DAF4A", 
                    "Breast" = "#984EA3")

treatment_colors <- c("Pre" = "#377EB8", "Post" = "#E41A1C")

print("Setup complete! All libraries loaded.")

