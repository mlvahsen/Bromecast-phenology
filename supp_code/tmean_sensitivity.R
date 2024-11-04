# Relative sensitivities analysis with common garden locations ONLY (i.e., no
# gravel treatments) compared to source climate for genotypes. Also creates
# Figure S7.

# Load libraries
library(tidyverse); library(here);
library(geosphere); library(usmap);
library(patchwork); library(ggnewscale);
library(here)

# Source in compiled data for the model
source(here("supp_code", "compile_data.R"))

# Set plotting theme
theme_set(theme_bw(base_size = 16))

# Read in climate of origin data for genotypes (tmean October - June). These
# data were downloaded via the prism_wrangle.R file in supp_code/
genotype_tmean <- read_csv("supp_data/genotype_tmean_norms.csv")

# Read in gps data for genotypes from main Bromecast repository
genotype_gps <- read_csv("https://raw.githubusercontent.com/pbadler/bromecast-data/main/gardens/rawdata/collection_sites.csv")
genotype_gps %>% 
  dplyr::select(Site.code = `Site code`, 
                longitude = Longitude,
                latitude = Latitude) -> genotype_gps

# Read in genotype codes from main Bromecast repository
genotypes <- read_csv("https://raw.githubusercontent.com/pbadler/bromecast-data/main/gardens/rawdata/sitecode2genotypenumber.csv")

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

# Get predicted means by genotype from linear mixed model
brms_lin <- read_rds("outputs/phenology_nokin_final.rds")
as_draws_df(brms_lin) -> obj

# Get genotype random intercept draws only
obj %>% 
  dplyr::select(contains("r_genotype") & contains("Intercept")) -> r_intercepts

# Rename columns with genotype number only
colnames(r_intercepts) <- parse_number(colnames(r_intercepts))

# Reorder columns to be in genotype number order
r_intercepts[,order(parse_number(colnames(r_intercepts)))] -> r_intercepts

# Get scaled values of PC axes
phen_flower_kin %>%
  mutate(pc1_sc = scale(phen_flower_kin$pc1)[,1],
         pc2_sc = scale(phen_flower_kin$pc2)[,1]) %>% 
  select(genotype, pc1_sc, pc2_sc) %>% 
  distinct() %>% 
  mutate(genotype = as.numeric(as.character(genotype))) %>% 
  arrange(genotype) -> pcs

# Get predicted means and CIs for each genotype
pred_all <- matrix(NA, nrow = length(obj$b_Intercept), ncol = length(unique(phen_flower_kin$genotype)))

for(i in 1:length(unique(phen_flower_kin$genotype))){
  pred_all[,i] <- obj$b_Intercept + as.data.frame(r_intercepts)[,i] + pcs$pc1_sc[i] * obj$b_pc1_sc + pcs$pc2_sc[i] * obj$b_pc2_sc
}

tibble(genotype = as.factor(parse_number(colnames(r_intercepts))),
       pred_mean = colMeans(pred_all),
       pred_upper = apply(pred_all, 2, quantile, probs = 0.975),
       pred_lower = apply(pred_all, 2, quantile, probs = 0.025),
       variance = apply(pred_all, 2, var)) %>% 
  merge(phen_flower_kin %>% select(genotype, site_code) %>% distinct()) %>% 
  merge(genotypes_tmean %>% select(site_code, tmean_mean)) %>% 
  mutate(diff = (pred_upper - pred_lower)^2)-> genotypes_tmean

# Read in tmeans for each site (file created in prism_wrangle.R)
site_tmean <- read_csv("supp_data/site_tmean.csv")

# Calculate relationship between CI and variance
ci_mod <- lm(variance ~ diff, data = genotypes_tmean)

# Get predicted means for each site gravel combo
site_preds <- sjPlot::plot_model(brms_lin, type = "emm", terms = "site")

# Format and merge with site mean temperature data (over the course of the
# growing season)
tibble(site_code = c("SS", "BA", "WI", "CH"),
       pred_mean = site_preds$data$predicted,
       pred_lower = site_preds$data$conf.low,
       pred_upper = site_preds$data$conf.high) %>%
  mutate(variance = coef(ci_mod)[1] + coef(ci_mod)[2] * (pred_upper - pred_lower)^2) %>% 
  merge(site_tmean) -> site_tmean_all

# Formatting genotypes data to merge
genotypes_tmean %>% 
  mutate(id = paste("genotype", genotype, sep = "_"),
         type = "Source environment") %>% 
  dplyr::select(jday = pred_mean, lower = pred_lower, upper = pred_upper,
                tmean = tmean_mean, id, type, variance) -> genotypes_tmean_tomerge

# Formatting sites data to merge
site_tmean_all %>% 
  mutate(id = site_code,
         type = "Current environment") %>% 
  dplyr::select(jday = pred_mean, lower = pred_lower, upper = pred_upper,
                tmean = mean_temp, id, type, variance) -> site_tmean_tomerge

# Bring together genotypes and sites data
climate_sens <- rbind(genotypes_tmean_tomerge, site_tmean_tomerge)

# Fit linear model with treatment contrasts and get CIs for interaction
options(contrasts = c("contr.treatment", "contr.poly"))
climate_sens$type <- factor(climate_sens$type, levels = c("Current environment", "Source environment"))
mod <- lm(jday ~ tmean * type, data = climate_sens, weights = 1/variance)
summary(mod) # Get coefficients
confint(mod) # Get confidence intervals

# Create Figure S7
text_x2 <- grid::textGrob("Avg. mean temperature October - June 1981-2010 (°C)",
                          gp=grid::gpar(fontsize=14, fontface="bold", col = "#2c7bb6"))

png("figs/FigS7_relative_sensitivies.png", width = 6, height = 5, res = 300, units = "in")
climate_sens %>% 
  ggplot(aes(x = tmean, y = jday, color = type, fill = type)) +
  geom_segment(aes(x = tmean, xend = tmean, y = lower, yend = upper), color = "black") +
  geom_point(size = 1) +
  geom_smooth(method = "lm", mapping = aes(weight = 1/variance)) + 
  geom_point(size = 1, shape = 21, color = "black", alpha = 0.8) +
  scale_color_manual(values = c("#fdae61", "#2c7bb6")) +
  scale_fill_manual(values = c("#fdae61", "#2c7bb6")) +
  labs(x = "Mean temperature October 2021 - June 2022 (°C)",
       y = "First day of flowering",
       color = "",
       fill = "") +
  theme(legend.position = "top",
        axis.title.x.bottom = element_text(color = "#fdae61", face = "bold", size = 14)) +
  theme(plot.margin = unit(c(1,1,2,1), "lines")) +
  annotation_custom(text_x2, xmin=0,xmax=16,ymin=80,ymax=80)+
  coord_cartesian(clip = "off") +
  scale_x_continuous(breaks = seq(3,12,by=3))+
  annotate("text", label = expression(paste(beta[temp:env], " = 5.8 (2.9, 8.7)")), x = 11, y = 205, size = 5) 
dev.off()
