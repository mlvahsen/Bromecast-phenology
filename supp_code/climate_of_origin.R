# This code calculates principal components analysis on bioclimate data for each
# genotype to get PC axes that describe the climate of origin for genotypes from
# the common garden site to use as a covariate in the phenology models.

## Preliminaries ####
# Load libraries
library(tidyverse); library(here);
library(geosphere); library(usmap);
library(patchwork); library(ggnewscale)

# Set plotting theme
theme_set(theme_bw(base_size = 14))

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
             size = 4, color = "dodgerblue", fill = "black", stroke = 1) +
  theme_classic(base_size = 14) +
  scale_fill_distiller(palette = "PiYG", limits = c(-6.8,6.8)) +
  labs(fill = "PC 1",
       x = "longitude",
       y = "latitude",
       shape = "common garden") +
  scale_shape_manual(values = c(0,5,2,6)) +
  ggtitle("(a) PC 1: cool & wet → hot & dry") -> pc1

# Plot of PC2 by site for genotypes
ggplot(data=region, aes(x=long, y=lat, group = group)) +
  geom_polygon(color = "gray", fill = "white")+
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(),
        axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
  coord_cartesian(ylim = c(34, 51), xlim = c(-125, -102)) +
  geom_point(data = genotypes_pc, aes(x = lon, y = lat, group = NA, fill = PC2),
             color = "black", shape = 21, size = 3) +
  geom_point(data = cg_pc, aes(x = lon, y = lat, group = NA, shape = site_code),
             size = 4, color = "dodgerblue", fill = "black", stroke = 1) +
  theme_classic(base_size = 14) +
  scale_fill_distiller(palette = "PuOr", limits = c(-4.75, 4.15)) +
  labs(fill = "PC 2",
       x = "longitude",
       y = "latitude",
       shape = "common garden") +
  scale_shape_manual(values = c(0,5,2,6)) +
  ggtitle("(b) PC 2: low → high seasonality") -> pc2

# Plot of common garden sites
cg_pc %>% 
  ggplot(aes(x = lon, y = elevation, fill = PC1, shape = site_code)) + 
  geom_point(size = 15) + 
  scale_fill_distiller(palette = "PiYG", limits = c(-6.8,6.8)) +
  guides(fill = "none") +
  new_scale_fill() +
  geom_point(aes(fill = PC2, shape = site_code), size = 6) + 
  scale_fill_distiller(palette = "PuOr", limits = c(-4.75,4.15)) +
  ylim(600, 2100) + xlim(-120, -102) +
  scale_shape_manual(values = 22:25) +
  ylab("elevation (m)") +
  xlab("longitude") +
  guides(shape = "none", fill = "none") +
  annotate(geom = "text", x = -117, y = 980, label = "Wildcat (WI)", size = 5) +
  annotate(geom = "text", x = -117, y = 1430, label = "Baltzor (BA)", size = 5) +
  annotate(geom = "text", x = -112, y = 1870, label = "Sheep Station (SS)", size = 5) +
  annotate(geom = "text", x = -105, y = 2090, label = "Cheyenne (CH)", size = 5) +
  geom_curve(aes(x = -110.5, y = 1480, xend = -112.2, yend = 1680), linewidth = 0.3,
             color = "gray37", arrow = arrow(length = unit(0.08, "inches")), curvature = -0.3) +
  geom_curve(aes(x = -108, y = 1620, xend = -111.6, yend = 1670), linewidth = 0.3, curvature = 0.3,
             color = "gray37", arrow = arrow(length = unit(0.08, "inches"))) +
  geom_text(aes(x = -108, y = 1480), label = "PC 2 value", fontface = "italic", color = "gray37") +
  geom_text(aes(x = -105.8, y = 1610), label = "PC 1 value", fontface = "italic", color = "gray37") +
  ggtitle("(c) common gardens") -> cg_pc_plot

png(here("figs/Fig1_ClimateOrigin.png"), height = 6.8, width = 13.3, res = 300, units = "in")
pc1 + pc2 +cg_pc_plot + plot_layout(guides = "collect", nrow = 1) & theme(legend.direction = "vertical",
                                                                legend.box = "horizontal", legend.position = "bottom")

dev.off()

## Create data frame for analysis ####

# Pull PC data for the first 2 axes and line up with genotype data
genotypes_pc %>% 
  dplyr::select(site_code, genotype, pc1 = PC1, pc2 = PC2) %>% 
  arrange(genotype) -> genotype_PCclimate

# Remove all objects except climate data
rm(list=setdiff(ls(), "genotype_PCclimate"))
