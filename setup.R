if (!dir.exists("figures")) dir.create("figures")
if (!dir.exists("tuning")) dir.create("tuning")

# Install pacman if it's not already installed
if (!require(pacman)) install.packages("pacman")

# Use pacman to install/load the required libraries for the SDM scripts
pacman::p_load(
  NicheMapR, magrittr, dplyr, tidyr, janitor, raster, sp, ggplot2, rgbif, dismo,
  ENMeval, sf, readr, terra, geodata, rnaturalearth, corrplot, gtools, igraph, usdm,
  factoextra, gridExtra, ecospat, predicts, readxl, stringr, RColorBrewer,
  tidyverse
)

library(SNPfiltR)
library(vcfR)
library(adegenet)
library(hierfstat)
library(pegas)
library(poppr)
library(PopGenReport)

# Change ggplot theme
theme_set(
  theme_classic() +
    theme(
      panel.border = element_rect(colour = "black", fill = NA),
      axis.text = element_text(colour = "black"),
      axis.title.x = element_text(margin = unit(c(2, 0, 0, 0), "mm")),
      axis.title.y = element_text(margin = unit(c(0, 4, 0, 0), "mm")),
      legend.position = "none"
    )
)

# Set the theme for the maps
theme_opts = list(
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    panel.background = element_rect(fill = 'white', colour = NA),
    plot.background = element_rect(),
    axis.line = element_blank(),
    axis.text.x = element_text(colour = "black"),
    axis.text.y = element_text(colour = "black"),
    axis.ticks = element_line(colour = "black"),
    axis.title.x = element_text(colour = "black"),
    axis.title.y = element_text(colour = "black"),
    plot.title = element_text(colour = "black"),
    panel.border = element_rect(fill = NA),
    legend.key = element_blank()
  )
)

