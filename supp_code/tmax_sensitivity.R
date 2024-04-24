# Relative sensitivities analysis with common garden locations AND gravel
# treatments compared to source climate for genotypes

# Load libraries
library(tidyverse); library(here);
library(geosphere); library(usmap);
library(patchwork); library(ggnewscale);
library(grid)

# Source in compiled data for the model
source(here("supp_code", "compile_data.R"))

# Set plotting theme
theme_set(theme_bw(base_size = 16))

# Read in climate of origin data for genotypes (tmax Jan - April)
genotype_tmax <- read_csv("supp_data/genotype_tmax_norms.csv")

# Get soil temperature data
temp <- read_csv("data/BCtemploggers.csv")

temp %>% 
  filter(Date > "2022-01-01" & Date < "2022-05-01") %>% 
  filter(Site != "Sheep" | Color != "White") -> clean_dat

soilT_to_airT0 <- function(x){
  pred <- 0.9164381*x + 8.838341
  return(pred)
}

soilT_to_airT0(10)

cbind(clean_dat, tmax = soilT_to_airT0(clean_dat$Temp_C)) %>% 
  group_by(Site, Color) %>% 
  summarize(tmax_mean = mean(tmax, na.rm = T)) %>% 
  ungroup() -> site_tmax

site_tmax %>% 
  mutate(site = case_when(Site == "Balzor" ~ "BA",
                          Site == "Cheyenne" ~ "CH",
                          Site == "Wildcat" ~ "WI",
                          Site == "Sheep" ~ "SS"),
         color = ifelse(Color == "Black", "black", "white")) %>% 
  dplyr::select(site, color, tmax_mean) -> site_tmax

# Read in garden gps data
garden_gps <- read_csv(here("~/Git/Bromecast/gardens/rawdata/garden_info.csv"))

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
  merge(genotype_tmax, all.y = T) -> genotypes_tmax

# Get predicted means by genotype
brms_lin <- read_rds("~/Desktop/phenology_nokin_final.rds")
as_draws_df(brms_lin) -> obj

obj %>% 
  select(contains("r_genotype") & contains("Intercept")) %>% 
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
  merge(genotypes_tmax) -> genotypes_tmax

# Get predicted means for each site gravel combo
site_gravel_preds <- sjPlot::plot_model(brms_lin, type = "emm", terms = c("site", "gravel"))

tibble(site = rep(c("SS", "BA", "WI", "CH"), each = 2),
       color = site_gravel_preds$data$group,
       jday = site_gravel_preds$data$predicted) %>% 
  merge(site_tmax) -> site_tmax_all

genotypes_tmax %>% 
  mutate(id = paste("genotype", genotype, sep = "_"),
         type = "Source environment") %>% 
  dplyr::select(jday = jday, tmax = tmax_mean, id, type) -> genotypes_tmax_tomerge

site_tmax_all %>% 
  mutate(id = paste(site, color, sep = "_"),
         type = "Current environment") %>% 
  dplyr::select(jday = jday, tmax = tmax_mean, id, type) -> site_tmax_tomerge

climate_sens <- rbind(genotypes_tmax_tomerge, site_tmax_tomerge)

# Fit linear model with treatment contrasts and get CIs for interaction
options(contrasts = c("contr.treatment", "contr.poly"))
climate_sens$type <- factor(climate_sens$type, levels = c("Current environment", "Source environment"))
mod <- lm(jday ~ tmax * type, data = climate_sens)
summary(mod)
confint(mod)

# Slope of local adaptation is the slope coefficient plus interaction term
coef(mod)[2] + coef(mod)[4]
# -1.036538
# One degree increase in temp means ~1 additional day earlier

# Slope of plasticity is the slope coefficient 
coef(mod)[2] 
# -5.374806
# One degree increase in temp means an additional 5 days earlier

text_x2 <- grid::textGrob("Avg. maximum temperature Jan - April 1981-2010 (°C)",
                          gp=gpar(fontsize=14, fontface="bold", col = "#2c7bb6"))

png("figs/Fig6_relative_sensitivies.png", width = 6, height = 5, res = 300, units = "in")
climate_sens %>% 
  ggplot(aes(x = tmax, y = jday, color = type, fill = type)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm") + 
  geom_point(size = 3, shape = 21, color = "black", alpha = 0.8) +
  scale_color_manual(values = c("#fdae61", "#2c7bb6")) +
  scale_fill_manual(values = c("#fdae61", "#2c7bb6")) +
  labs(x = "Maximum temperature Jan - April 2022 (°C)",
       y = "First day of flowering",
       color = "",
       fill = "") +
  theme(legend.position = "top",
        axis.title.x.bottom = element_text(color = "#fdae61", face = "bold", size = 14)) +
  theme(plot.margin = unit(c(1,1,2,1), "lines")) +
  annotation_custom(text_x2, xmin=0,xmax=21.5,ymin=103,ymax=103)+
  ylim(118, 175) +
  coord_cartesian(clip = "off") +
  annotate("text", label = expression(paste(beta[temp:env], " = 4.3 (2.7, 6.0)")), x = 15.2, y = 175, size = 5) 
dev.off()
