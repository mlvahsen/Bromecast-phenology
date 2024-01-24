# Relative sensitivities analysis with common garden locations ONLY (i.e., no
# gravel treatments) compared to source climate for genotypes

# Load libraries
library(tidyverse); library(here);
library(geosphere); library(usmap);
library(patchwork); library(ggnewscale);
library(here)

# Source in compiled data for the model
source(here("supp_code", "compile_data.R"))

# Set plotting theme
theme_set(theme_bw(base_size = 16))

# Read in climate of origin data for genotypes (tmean October - June)
genotype_tmean <- read_csv("~/Desktop/genotype_tmean_norms.csv")

# Read in gps data for genotypes
genotype_gps <- read_csv(here("~/Git/Bromecast/gardens/rawdata/collection_sites.csv"))
genotype_gps %>% 
  dplyr::select(Site.code = `Site code`, 
                longitude = Longitude,
                latitude = Latitude) -> genotype_gps

# Read in genotype codes
genotypes <- read_csv(here("~/Git/Bromecast/gardens/rawdata/sitecode2genotypenumber.csv"))

# Merge together genotype gps and code data and have just one row per genotype
merge(genotype_gps, genotypes) %>% 
  dplyr::select(site_code = Site.code,
                longitude,
                latitude,
                genotypeID) %>% 
  distinct() %>% 
  mutate(genotype = parse_number(genotypeID)) %>% 
  arrange(genotype) -> genotypes_with_gps

# Add genotype number to the data
genotypes_with_gps %>% 
  dplyr::select(genotype, site_code) %>% 
  merge(genotype_tmean, all.y = T) -> genotypes_tmean

# Get predicted means by genotype
brms_lin <- read_rds("~/Desktop/brms_linear.rds")
as_draws_df(brms_lin) -> obj

obj %>% 
  dplyr::select(contains("r_genotype") & contains("Intercept")) %>% 
  gather(key = genotype, value = intercept) %>% 
  mutate(genotype = parse_number(genotype)) %>% 
  group_by(genotype) %>% 
  summarize(mean = mean(intercept)) %>% 
  arrange(genotype) %>% ungroup()-> random_intercepts

phen_flower_kin %>%
  mutate(pc1_sc = scale(phen_flower_kin$pc1)[,1],
         pc2_sc = scale(phen_flower_kin$pc2)[,1]) %>% 
  group_by(genotype) %>% 
  summarize(pc1_sc = mean(pc1_sc),
            pc2_sc = mean(pc2_sc)) %>% 
  mutate(genotype = as.numeric(as.character(genotype))) %>% 
  arrange(genotype) %>% ungroup() -> pcs

pred <- mean(obj$b_Intercept) + random_intercepts$mean + pcs$pc1_sc * mean(obj$b_pc1_sc) + pcs$pc2_sc * mean(obj$b_pc2_sc)

# Store parameters to unscale later
mean_pc1 <- attr(scale(phen_flower_kin$pc1),"scaled:center")
sd_pc1 <- attr(scale(phen_flower_kin$pc1),"scaled:scale")
mean_pc2 <- attr(scale(phen_flower_kin$pc2),"scaled:center")
sd_pc2 <- attr(scale(phen_flower_kin$pc2),"scaled:scale")

cbind(pcs, jday = pred) -> preds_by_genotype
preds_by_genotype$pc1 <- preds_by_genotype$pc1_sc * sd_pc1 + mean_pc1
preds_by_genotype$pc2 <- preds_by_genotype$pc2_sc * sd_pc2 + mean_pc2

preds_by_genotype %>% 
  merge(genotypes_tmean) -> genotypes_tmean

# Get predicted means for each site gravel combo
site_preds <- sjPlot::plot_model(brms_lin, type = "emm", terms = c("site"))

# Read in tmeans for each site
site_tmean <- read_csv("~/Desktop/site_tmean.csv")

tibble(site_code = rep(c("SS", "BA", "WI", "CH")),
       jday = site_preds$data$predicted) %>% 
  merge(site_tmean) -> site_tmean_all

genotypes_tmean %>% 
  mutate(id = paste("genotype", genotype, sep = "_"),
         type = "Local adaptation") %>% 
  dplyr::select(jday = jday, tmean = tmean_mean, id, type) -> genotypes_tmean_tomerge

site_tmean_all %>% 
  mutate(id = site_code,
         type = "Plasticity") %>% 
  dplyr::select(jday = jday, tmean = mean_temp, id, type) -> site_tmean_tomerge

climate_sens <- rbind(genotypes_tmean_tomerge, site_tmean_tomerge)

# Fit linear model with treatment contrasts and get CIs for interaction
options(contrasts = c("contr.treatment", "contr.poly"))
climate_sens$type <- factor(climate_sens$type, levels = c("Plasticity", "Local adaptation"))
mod <- lm(jday ~ tmean * type, data = climate_sens)
summary(mod)
confint(mod)

text_x2 <- grid::textGrob("Avg. mean temperature October - June 1981-2010 (°C)",
                          gp=grid::gpar(fontsize=14, fontface="bold", col = "#2c7bb6"))

png("figs/FigS7_relative_sensitivies.png", width = 6, height = 5, res = 300, units = "in")
climate_sens %>% 
  ggplot(aes(x = tmean, y = jday, color = type, fill = type)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm") + 
  geom_point(size = 3, shape = 21, color = "black", alpha = 0.8) +
  scale_color_manual(values = c("#fdae61", "#2c7bb6")) +
  scale_fill_manual(values = c("#fdae61", "#2c7bb6")) +
  labs(x = "Mean temperature October 2021 - June 2022 (°C)",
       y = "First day of flowering",
       color = "",
       fill = "") +
  theme(legend.position = "top",
        axis.title.x.bottom = element_text(color = "#fdae61", face = "bold", size = 14)) +
  theme(plot.margin = unit(c(1,1,2,1), "lines")) +
  annotation_custom(text_x2, xmin=0,xmax=17.5,ymin=77,ymax=77)+
  coord_cartesian(clip = "off") +
  annotate("text", label = expression(paste(beta[temp:type], " = 5.9 (1.7, 10.0)")), x = 11.9, y = 205, size = 5) 
dev.off()
