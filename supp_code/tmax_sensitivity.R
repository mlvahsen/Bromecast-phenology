# Relative sensitivities analysis with common garden locations AND gravel
# treatments compared to source climate for genotypes. Also makes Figure 6.

# Load libraries
library(tidyverse); library(here);
library(geosphere); library(usmap);
library(patchwork); library(ggnewscale);
library(grid)

# Source in compiled data for the model
source(here("supp_code", "compile_data.R"))

# Set plotting theme
theme_set(theme_bw(base_size = 16))

# Read in climate of origin data for genotypes (tmax October - June). This is
# created from the prism_wrangle.R script.
genotype_tmax <- read_csv("supp_data/genotype_tmax_norms.csv")

# Get soil temperature data
temp <- read_csv("data/BCtemploggers.csv")

# Restrict to Jan - Apr 2021
temp %>% 
  filter(Date > "2022-01-01" & Date < "2022-05-01") -> clean_dat

# Create soil to air temperature function based on supplemental figure from
# Petrie et al. (2020, Climate Change)
soilT_to_airT0 <- function(x){
  pred <- 0.9164381*x + 8.838341
  return(pred)
}

# Calculate cg site and gravel air temperatures 
cbind(clean_dat, tmax = soilT_to_airT0(clean_dat$Temp_C)) %>% 
  group_by(Site, Color) %>% 
  summarize(tmax_mean = mean(tmax, na.rm = T)) %>% 
  ungroup() -> site_tmax

# Tidy up data
site_tmax %>% 
  mutate(site = case_when(Site == "Balzor" ~ "BA",
                          Site == "Cheyenne" ~ "CH",
                          Site == "Wildcat" ~ "WI",
                          Site == "Sheep" ~ "SS"),
         color = ifelse(Color == "Black", "black", "white")) %>% 
  dplyr::select(site, color, tmax_mean) -> site_tmax

# Read in garden gps data from main Bromecast repository
garden_gps <- read_csv("https://raw.githubusercontent.com/pbadler/bromecast-data/main/gardens/rawdata/garden_info.csv")

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
  merge(genotype_tmax, all.y = T) -> genotypes_tmax

# Get predicted means by genotype - use model output from LMM
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
  merge(genotype_tmax %>% select(site_code, tmax_mean)) -> genotype_tmax

# Calculate relationship between CI and variance
ci_mod <- lm(variance ~ diff, data = genotype_tmax %>% mutate(diff = (pred_upper - pred_lower)^2))

# Get predicted means for each site gravel combo
site_gravel_preds <- sjPlot::plot_model(brms_lin, type = "emm", terms = c("site", "gravel"))

# Format and merge with site max temperature data (over the course of the
# growing season)
tibble(site = rep(c("SS", "BA", "WI", "CH"), each = 2),
       color = site_gravel_preds$data$group,
       pred_mean = site_gravel_preds$data$predicted,
       pred_lower = site_gravel_preds$data$conf.low,
       pred_upper = site_gravel_preds$data$conf.high) %>%
  mutate(variance = coef(ci_mod)[1] + coef(ci_mod)[2] * (pred_upper - pred_lower)^2) %>% 
  merge(site_tmax) -> site_tmax_all

# Clean up genotype data so it can be merged with cg site data
genotype_tmax %>% 
  mutate(id = paste("genotype", genotype, sep = "_"),
         type = "Source environment") %>% 
  dplyr::select(jday = pred_mean, tmax = tmax_mean, id, type, variance,
                upper = pred_upper, lower = pred_lower) -> genotypes_tmax_tomerge

# Clean up site data so it can be merged with genotype data
site_tmax_all %>% 
  mutate(id = paste(site, color, sep = "_"),
         type = "Current environment") %>% 
  dplyr::select(jday = pred_mean, tmax = tmax_mean, id, type, variance,
                upper = pred_upper, lower = pred_lower) -> site_tmax_tomerge

# Merge site and genotype info into one dataset
climate_sens <- rbind(genotypes_tmax_tomerge, site_tmax_tomerge)

# Add code for different sites and different gravel treatments
climate_sens %>% 
  mutate(site_code = case_when(grepl("BA", id) ~ "Cool, aseasonal (BA)",
                               grepl("SS", id) ~ "Cold, aseasonal (SS)",
                               grepl("WI", id) ~ "Hot, seasonal (WI)",
                               grepl("CH", id) ~ "Cool, seasonal (CH)",
                               T ~ "genotype")) %>% 
  mutate(gravel_code = case_when(grepl("black", id) ~ "High temp.\n(black gravel)",
                                 grepl("white", id) ~ "Low temp.\n(white gravel)",
                                 T ~ "genotype")) -> climate_sens

# Fit weighted linear model with treatment contrasts and get CIs for
# interaction. Weights are the inverse of variance.
options(contrasts = c("contr.treatment", "contr.poly"))
climate_sens$type <- factor(climate_sens$type, levels = c("Current environment", "Source environment"))
mod <- lm(jday ~ tmax * type, data = climate_sens, weights = 1/variance)
summary(mod) # Get coefficients
confint(mod) # Get confidence intervals

# Slope of local adaptation is the slope coefficient plus interaction term
coef(mod)[2] + coef(mod)[4]
# -1.030154 
# One degree increase in temp means ~1 additional day earlier

# Slope of plasticity is the slope coefficient 
coef(mod)[2] 
# -5.288871 
# One degree increase in temp means an additional ~5 days earlier

# Make Figure 6
text_x2 <- grid::textGrob("Avg. maximum temperature Jan - Apr 1981-2010 (°C)",
                          gp=gpar(fontsize=14, fontface="bold", col = "#2c7bb6"))

climate_sens %>% 
  ggplot(aes(x = tmax, y = jday, color = type, fill = type)) +
  geom_segment(aes(x = tmax, xend = tmax, y = lower, yend = upper), color = "black") +
  geom_point(size = 2) +
  geom_smooth(method = "lm", mapping = aes(weight = 1/variance)) + 
  geom_point(size = 2, shape = 21, color = "black", alpha = 0.8) +
  scale_color_manual(values = c("#fdae61", "#2c7bb6")) +
  scale_fill_manual(values = c("#fdae61", "#2c7bb6")) +
  labs(x = "Maximum temperature Jan 2022 - Apr 2022 (°C)",
       y = "First day of flowering",
       color = "",
       fill = "") +
  guides(shape = FALSE) +
  theme(legend.position = "top",
        axis.title.x.bottom = element_text(color = "#fdae61", face = "bold", size = 14)) +
  theme(plot.margin = unit(c(1,1,2,1), "lines")) +
  annotation_custom(text_x2, xmin=0,xmax=21.5,ymin=103,ymax=103)+
  ylim(118, 175) +
  coord_cartesian(clip = "off") +
  annotate("text", label = expression(paste(beta[temp:env], " = 4.3 (3.2, 5.4)")), x = 15.2, y = 175, size = 5) -> Fig6

ggsave("figs/Figure6.eps", Fig6, width = 6, height = 5, units = "in", dpi = 600, device = cairo_ps)

