# This code calculates principal components analysis on bioclimate data for each
# genotype to get PC axes that describe the climate of origin for genotypes from
# the common garden site to use as a covariate in the phenology models.

## Preliminaries ####
# Load libraries
library(tidyverse); library(here);
library(geosphere); library(usmap);
library(patchwork)

# Read in bioclimatic data for genotypes (and common garden locations)
bioclim <- read_csv("~/Documents/Git/Bromecast Data/gardens/deriveddata/BioclimateOfOrigin_AllGenotypes.csv")

## Run PCA ####

# First center and scale all bioclimatic variables
bioclim[,6:ncol(bioclim)] <- apply(bioclim[,6:ncol(bioclim)], 2, scale)

# Run PCA
pca_out <- prcomp(bioclim[,6:ncol(bioclim)])

# Get percent explained by each PC axis
round(pca_out$sdev^2 / sum(pca_out$sdev^2),3) -> perc_explained

# Bind PC axis data with original data
cbind(bioclim, pca_out$x) -> bioclim_pc

## Create maps of PC variables ####

# Get US map with states
state <- map_data("state")
# Get Canada map
canada <- map_data("world") %>% filter(region == "Canada")
# Join into single region
region <- rbind(state, canada)

# Divide dataset into genotype collection sites and common garden sites
`%notin%` <- Negate(`%in%`)
genotypes_pc <- bioclim_pc %>% filter(site_code %notin% c("SS", "CH", "BA", "WI"))
cg_pc <- bioclim_pc %>% filter(site_code %in% c("SS", "CH", "BA", "WI"))

# Plot of elevation by site for genotypes
ggplot(data=region, aes(x=long, y=lat, group = group)) +
  geom_polygon(color = "gray", fill = "white")+
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(),
        axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
  coord_cartesian(ylim = c(34, 51), xlim = c(-125, -102)) +
  geom_point(data = genotypes_pc, aes(x = lon, y = lat, group = NA, fill = elevation),
             color = "black", shape = 21, size = 3) +
  geom_point(data = cg_pc, aes(x = lon, y = lat, group = NA, shape = site_code),
             size = 2, color = "dodgerblue", fill = "black", stroke = 1) +
  theme_classic(base_size = 14) +
  scale_fill_distiller(palette = "BrBG") +
  labs(fill = "elevation (m)",
       x = "longitude",
       y = "latitude",
       shape = "common garden") +
  scale_shape_manual(values = 22:25) +
  ggtitle("(a) elevation")-> elevation

# Plot of PC1 by site for genotypes
ggplot(data=region, aes(x=long, y=lat, group = group)) +
  geom_polygon(color = "gray", fill = "white")+
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(),
        axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
  coord_cartesian(ylim = c(34, 51), xlim = c(-125, -102)) +
  geom_point(data = genotypes_pc, aes(x = lon, y = lat, group = NA, fill = PC1),
             color = "black", shape = 21, size = 3) +
  geom_point(data = cg_pc, aes(x = lon, y = lat, group = NA, shape = site_code),
             size = 2, color = "dodgerblue", fill = "black", stroke = 1) +
  theme_classic(base_size = 14) +
  scale_fill_distiller(palette = "PiYG") +
  labs(fill = "PC 1",
       x = "longitude",
       y = "latitude",
       shape = "common garden") +
  scale_shape_manual(values = 22:25) +
  ggtitle("(b) PC 1")-> pc1

# Plot of PC2 by site for genotypes
ggplot(data=region, aes(x=long, y=lat, group = group)) +
  geom_polygon(color = "gray", fill = "white")+
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(),
        axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
  coord_cartesian(ylim = c(34, 51), xlim = c(-125, -102)) +
  geom_point(data = genotypes_pc, aes(x = lon, y = lat, group = NA, fill = PC2),
             color = "black", shape = 21, size = 3) +
  geom_point(data = cg_pc, aes(x = lon, y = lat, group = NA, shape = site_code),
             size = 2, color = "dodgerblue", fill = "black", stroke = 1) +
  theme_classic(base_size = 14) +
  scale_fill_distiller(palette = "PuOr") +
  labs(fill = "PC 2",
       x = "longitude",
       y = "latitude",
       shape = "common garden") +
  scale_shape_manual(values = 22:25) +
  ggtitle("(c) PC 2") -> pc2

png(here("figs/FigSXX_ClimateOrigin.png"), height = 6, width = 10, res = 300, units = "in")
elevation + pc1 + pc2 + plot_layout(guides = "collect") & theme(legend.direction = "vertical",
                                                                legend.box = "horizontal", legend.position = "bottom")

dev.off()

## Create data frame for analysis ####

# Pull PC data for the first 2 axes and line up with genotype data
genotypes_pc %>% 
  dplyr::select(site_code, genotype, pc1 = PC1, pc2 = PC2) %>% 
  arrange(genotype) -> genotype_PCclimate

# Remove all objects except climate data
rm(list=setdiff(ls(), "genotype_PCclimate"))
