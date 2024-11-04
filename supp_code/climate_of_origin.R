# This code calculates principal components analysis on bioclimate data for each
# genotype to get PC axes that describe the climate of origin for genotypes from
# the common garden site to use as a covariate in the phenology models. Also
# creates Figure 2

## Preliminaries ####
# Load libraries
library(tidyverse); library(here);
library(geosphere); library(usmap);
library(patchwork); library(ggnewscale);
library(factoextra)

# Set plotting theme
theme_set(theme_bw(base_size = 14))

# Read in bioclimatic data for genotypes (and common garden locations) from main
# Bromecast repository
bioclim <- read_csv("https://raw.githubusercontent.com/pbadler/bromecast-data/main/gardens/deriveddata/BioclimateOfOrigin_AllGenotypes.csv")

## Run PCA ####

# Center and scale bioclimatic covariates
bioclim[,6:ncol(bioclim)] <- apply(bioclim[,6:ncol(bioclim)], 2, scale)

# Run PCA
pca_out <- prcomp(bioclim[,6:ncol(bioclim)])

# Get percent explained by each PC axis
round(pca_out$sdev^2 / sum(pca_out$sdev^2),3) -> perc_explained

# Get which variables are contributing most to each PC axis
get_pca_var(pca_out)$contrib[,1]
sort(get_pca_var(pca_out)$contrib[,2], decreasing = T)

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
cg_pc %>% 
  mutate(site_code = case_when(site_code == "SS" ~ "Cold, aseasonal (SS)",
                               site_code == "BA" ~ "Cool, aseasonal (BA)",
                               site_code == "CH" ~ "Cool, seasonal (CH)",
                               site_code == "WI" ~ "Hot, seasonal (WI)")) -> cg_pc

# Plot of PC1 by site for genotypes
ggplot(data=region, aes(x=long, y=lat, group = group)) +
  geom_polygon(color = "gray", fill = "white")+
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(),
        axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
  coord_cartesian(ylim = c(34, 51), xlim = c(-125, -102)) +
  geom_point(data = genotypes_pc, aes(x = lon, y = lat, group = NA, fill = PC1),
             color = "black", shape = 21, size = 3) +
  geom_point(data = cg_pc, aes(x = lon, y = lat, group = NA, shape = site_code),
             size = 4, color = "black", stroke = 1) +
  theme_classic(base_size = 18) +
  theme_classic(base_size = 14) +
  scale_fill_distiller(palette = "PiYG", limits = c(-7,7)) +
  labs(fill = "PC 1",
       x = "Longitude",
       y = "Latitude",
       shape = "Common garden") +
  scale_shape_manual(values = c(2,5,0,6)) +
  ggtitle(paste("(a) PC 1: cool & wet", "\U2192", "hot & dry")) -> pc1

# Plot of PC2 by site for genotypes
ggplot(data=region, aes(x=long, y=lat, group = group)) +
  geom_polygon(color = "gray", fill = "white")+
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(),
        axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
  coord_cartesian(ylim = c(34, 51), xlim = c(-125, -102)) +
  geom_point(data = genotypes_pc, aes(x = lon, y = lat, group = NA, fill = PC2),
             color = "black", shape = 21, size = 3) +
  geom_point(data = cg_pc, aes(x = lon, y = lat, group = NA, shape = site_code),
             stroke = 1, color = "black", size = 4) +
  theme_classic(base_size = 18) +
  theme_classic(base_size = 14) +
  scale_fill_distiller(palette = "PuOr", limits = c(-5, 5)) +
  labs(fill = "PC 2",
       x = "Longitude",
       y = "Latitude",
       shape = "Common garden") +
  scale_shape_manual(values = c(2,5,0,6)) +
  ggtitle(paste("(b) PC 2: low", "\U2192", "high seasonality")) -> pc2

# Plot of common garden sites
cg_pc %>%
  ggplot(aes(x = lon, y = elevation)) +
  geom_point(size = 15, aes(fill = PC1, shape = site_code)) +
  scale_fill_distiller(palette = "PiYG", limits = c(-7,7)) +
  guides(fill = "none") +
  new_scale_fill() +
  geom_point(aes(fill = PC2, shape = site_code), size = 6) +
  scale_fill_distiller(palette = "PuOr", limits = c(-5,5)) +
  ylim(120, 2406) + xlim(-121, -101) +
  scale_shape_manual(values = c(24,23,22,25)) +
  ylab("Elevation (m)") +
  xlab("Longitude") +
  guides(shape = "none", fill = "none") +
  geom_point(data = genotypes_pc, aes(x = lon, y = elevation), alpha = 0.2) +
  annotate(geom = "text", x = -117, y = 1030, label = "Hot, seasonal (WI)", size = 4) +
  annotate(geom = "text", x = -117, y = 1480, label = "Cool, aseasonal (BA)", size = 4) +
  annotate(geom = "text", x = -112, y = 1960, label = "Cold, aseasonal (SS)", size = 4) +
  annotate(geom = "text", x = -105, y = 2140, label = "Cool, seasonal (CH)", size = 4) +
  geom_curve(aes(x = -110.5, y = 1480, xend = -112.2, yend = 1680), linewidth = 0.3,
             color = "gray37", arrow = arrow(length = unit(0.08, "inches")), curvature = -0.3) +
  geom_curve(aes(x = -108, y = 1620, xend = -111.4, yend = 1670), linewidth = 0.3, curvature = 0.3,
             color = "gray37", arrow = arrow(length = unit(0.08, "inches"))) +
  geom_text(aes(x = -108, y = 1480), label = "PC 2 value", fontface = "italic", color = "gray37") +
  geom_text(aes(x = -105.8, y = 1610), label = "PC 1 value", fontface = "italic", color = "gray37") +
  ggtitle("(c) Common garden sites") -> cg_pc_plot

Fig2 <- pc1 + pc2 +cg_pc_plot + plot_layout(guides = "collect", nrow = 1) & theme(legend.direction = "vertical",
                                                                legend.box = "horizontal", legend.position = "bottom")

ggsave("figs/Figure2.eps", Fig2, width = 13.3, height = 6.8, units = "in", dpi = 600, device = cairo_ps)

## Create data frame for analysis ####

# Pull PC data for the first 2 axes and line up with genotype data
genotypes_pc %>% 
  dplyr::select(site_code, genotype, pc1 = PC1, pc2 = PC2) %>% 
  arrange(genotype) -> genotype_PCclimate

# Remove all objects except climate data -- this code gets sourced for main
# analysis
rm(list=setdiff(ls(), "genotype_PCclimate"))
