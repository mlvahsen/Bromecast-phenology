## Calculate with PC1 info ####

# Genotype ####
predicted_means %>% 
  group_by(pc1, genotype) %>% 
  summarize(mean_pred = mean(jday_pred)) %>% 
  ungroup() -> genotype_level_means

genotype_level_means %>% 
  ggplot(aes(x = pc1, y = log(mean_pred))) +
  geom_point(size = 4, shape = 21, fill = "gray") +
  geom_smooth(method = "lm")

genotype_mod <- lm(mean_pred ~ pc1, data = genotype_level_means)

# Common garden site ####

# Read in bioclimatic data for genotypes (and common garden locations)
bioclim <- read_csv("~/Git/Bromecast/gardens/deriveddata/BioclimateOfOrigin_AllGenotypes.csv")

# Run PCA 

# First center and scale all bioclimatic variables
bioclim[,6:ncol(bioclim)] <- apply(bioclim[,6:ncol(bioclim)], 2, scale)

# Run PCA
pca_out <- prcomp(bioclim[,6:ncol(bioclim)])

# Get percent explained by each PC axis
round(pca_out$sdev^2 / sum(pca_out$sdev^2),3) -> perc_explained

# Bind PC axis data with original data
cbind(bioclim, pca_out$x) -> bioclim_pc

# Divide dataset to common garden sites
cg_pc <- bioclim_pc %>% filter(site_code %in% c("SS", "CH", "BA", "WI"))

predicted_means %>% 
  group_by(site) %>% 
  summarize(mean_pred = mean(jday_pred)) %>% 
  ungroup() %>% 
  mutate(pc1 = cg_pc$PC1) -> site_level_means

site_level_means %>% 
  ggplot(aes(x = pc1, y = log(mean_pred))) +
  geom_point(size = 4, shape = 21, fill = "gray") +
  geom_smooth(method = "lm")

site_mod <- lm(mean_pred ~ pc1, data = site_level_means)

# Make joined plot ####
site_level_means %>% 
  select(index = site, mean_pred, pc1) -> site_level_means
genotype_level_means %>% 
  select(index = genotype, mean_pred, pc1) -> genotype_level_means

rbind(site_level_means, genotype_level_means) %>% 
  mutate(type = ifelse(index %in% c("BA", "CH", "SS", "WI"),
                       "Current climate", "Climate of origin")) -> joint_data

png("figs/Fig5_relative_sensitivies.png", width = 5.5, height = 5, res = 300, units = "in")
joint_data %>% 
  ggplot(aes(x = pc1, y = mean_pred + 108, color = type, fill = type)) +
  geom_smooth(method = "lm") + 
  geom_point(size = 3, shape = 21, color = "black", alpha = 0.8) +
  scale_color_manual(values = c("#2c7bb6", "#fdae61")) +
  scale_fill_manual(values = c("#2c7bb6", "#fdae61")) +
  xlab("PC 1 (cool & wet → hot & dry)") +
  ylab("Predicted first day of flowering") +
  labs(color = "Type", fill = "Type") +
  annotate("text", label = expression(paste(beta[PC1:type], " = (-5.41, -0.60)")), x = 4, y = 195, size = 5) +
  theme(legend.position = "top")
dev.off()

options(contrasts = c("contr.treatment", "contr.sum"))
joint_mod <- lm(mean_pred ~ pc1 * type, data = joint_data)
confint(joint_mod)

## Calculate with tmean info from PRISM data ####
site_temps <- tibble(site = c("BA", "CH", "SS", "WI"),
                     tmean = c(5.25, 5.33, 4.34, 6.73))

genotype_temps <- read_csv("~/Desktop/genotype_tmean_norms.csv")
genotypes <- read_csv("~/Git/Bromecast/gardens/rawdata/sitecode2genotypenumber.csv")
genotypes %>% 
  mutate(genotype = parse_number(genotypeID),
         site_code = Site.code) %>% 
  dplyr::select(genotype, site_code) %>% 
  merge(genotype_temps) -> genotype_temps

predicted_means %>% 
  group_by(genotype) %>% 
  summarize(mean = mean(jday_pred)+108) %>% 
  ungroup() %>% 
  merge(genotype_temps) -> genotype_temps_pred

predicted_means %>% 
  group_by(site) %>% 
  summarize(mean = mean(jday_pred)+108) %>% 
  ungroup() %>% 
  merge(site_temps) -> site_temps_pred

genotype_temps_pred %>% 
  dplyr::select(id = genotype,
         jday = mean,
         tmean = tmean_mean) %>% 
  mutate(type = "Source climate") -> genotype_temps_pred

site_temps_pred %>% 
  dplyr::select(id = site,
                jday = mean,
                tmean) %>% 
  mutate(type = "Current climate") -> site_temps_pred

rbind(genotype_temps_pred, site_temps_pred) -> joint_data

options(contrasts = c("contr.treatment", "contr.sum"))
joint_mod <- lm(jday ~ tmean * type, data = joint_data)
confint(joint_mod)

png("figs/Fig5_relative_sensitivies_tmean.png", height = 5, width = 6, res = 300, units = "in")
joint_data %>% 
  ggplot(aes(x = tmean, y = jday, color = type, fill = type)) +
  geom_smooth(method = "lm") + 
  geom_point(size = 3, shape = 21, color = "black", alpha = 0.8) +
  scale_color_manual(values = c("#fdae61", "#2c7bb6")) +
  scale_fill_manual(values = c("#fdae61", "#2c7bb6")) +
  xlab("Average mean temperature (°C) Oct - June") +
  ylab("Predicted first day of flowering") +
  labs(color = "Type", fill = "Type") +
  annotate("text", label = expression(paste(beta[temp:type], " = 10.0 (4.6, 15.4)")), x = 12, y = 200, size = 5) +
  theme(legend.position = "top") +
  ylim(105,200)
dev.off()
