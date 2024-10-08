# This code fits the linear models for flowering time: both version with and
# without kinship matrix included. This code also creates Figs 3, 4, 5, S4, S5

# Load libraries
library(tidyverse); library(mgcv); library(gratia); library(geomtextpath);
library(here); library(readr); library(brms); library(RcppCNPy); library(lme4);
library(patchwork); library(bayesplot); library(usmap); library(egg);
library(RVAideMemoire)

# Source in compiled data for the model
source(here("supp_code", "compile_data.R"))

phen_flower_kin %>% filter(is.na(resurrection_date)) -> phen_flower_kin

# Center and scale continuous variables
phen_flower_kin$pc1_sc <- scale(phen_flower_kin$pc1)[,1]
phen_flower_kin$pc2_sc <- scale(phen_flower_kin$pc2)[,1]

# Store parameters to unscale later
mean_pc1 <- attr(scale(phen_flower_kin$pc1),"scaled:center")
sd_pc1 <- attr(scale(phen_flower_kin$pc1),"scaled:scale")
mean_pc2 <- attr(scale(phen_flower_kin$pc2),"scaled:center")
sd_pc2 <- attr(scale(phen_flower_kin$pc2),"scaled:scale")

## Fit Bayesian model ####
# Fit Bayesian linear model

start <- Sys.time()

# Fit with summed contrasts to more easily infer global means
options(contrasts = c("contr.sum", "contr.poly"))

# Fit version of the model with kinship matrix
brms_lin <- brm(
  jday ~ 0 + Intercept + density * gravel * pc1_sc + density * gravel * pc2_sc +
    site * pc1_sc + site * pc2_sc + (1 + density + gravel + site || gr(genotype, cov = Amat)) +
    (1 | plot_unique),
  data = phen_flower_kin,
  data2 = list(Amat = kin),
  family = gaussian(),
  chains = 3, cores = 1, iter = 7500,
  # Set seed for reproducibility
  seed = 4685
)

# Fit version of the model with out kinship matrix
brms_lin_nokin <- brm(
  jday ~ 0 + Intercept + density * gravel * pc1_sc + density * gravel * pc2_sc +
    site * pc1_sc + site * pc2_sc + (1 + density + gravel + site || genotype) +
    (1 | plot_unique),
  data = phen_flower_kin,
  family = gaussian(),
  chains = 3, cores = 1, iter = 7500,
  # Set seed for reproducibility
  seed = 4685
)

# Fit version of the model with out kinship matrix and without source climate to
# quantify how much variance it explains
brms_lin_nokin_nosource <- brm(
  jday ~ 0 + Intercept + density * gravel + site +
    (1 + density + gravel + site || genotype) +
    (1 | plot_unique),
  data = phen_flower_kin,
  family = gaussian(),
  chains = 3, cores = 1, iter = 2000,
  # Set seed for reproducibility
  seed = 4685
)

# Fit version of the model with out kinship matrix and without current climate to
# quantify how much variance it explains
brms_lin_nokin_nocurrent <- brm(
  jday ~ 0 + Intercept + pc1 + pc2 + (1 | genotype) + (1 | plot_unique),
  data = phen_flower_kin,
  data2 = list(Amat = kin),
  family = gaussian(),
  chains = 3, cores = 1, iter = 2000,
  # Set seed for reproducibility
  seed = 4685
)
# Save model objects for use later
write_rds(brms_lin, "outputs/phenology_kin_final.rds")
write_rds(brms_lin_nokin, "outputs/phenology_nokin_final.rds")

# Read in Rdata object of no kinship matrix model (model that we are doing most
# of our inference from)
brms_lin_nokin <- read_rds("outputs/phenology_nokin_final.rds")

## Create graphics - Model checking ####

# Generate posterior predictive distribution
ppreds <- posterior_predict(brms_lin_nokin, draws = 1000)
# Check against distribution of data, predicted mean, and predicted SD
ppc_dens_overlay(phen_flower_kin$jday, ppreds)
# Fits data distribution reasonably well (smooths over high density sampling
# dates which is to be expected)
ppc_stat(phen_flower_kin$jday, ppreds, stat = "mean") -> ppc_mean
# Captures mean well
ppc_stat(phen_flower_kin$jday, ppreds, stat = "sd") -> ppc_sd

ppc_mean + xlab("Mean of posterior predictive distribution") +
  theme(legend.position = "none")->ppc_mean
ppc_sd + xlab("Std. dev. of posterior predictive distribution") +
  theme(legend.position = "none")->ppc_sd
# Capture standard deviation well

# Make a plot of predicted and observed for model
# Calculate R2 of observed vs predicted
give_me_R2 <- function(preds,actual){
  rss <- sum(( preds - actual ) ^ 2)  ## residual sum of squares
  tss <- sum((actual - mean(actual)) ^ 2)  ## total sum of squares
  rsq <- 1 - rss/tss
  return(rsq)
}

r2 <- give_me_R2(colMeans(ppreds), phen_flower_kin$jday)

# Create observed vs predicted plot and add R2 value
tibble(predicted = colMeans(ppreds),
       observed = phen_flower_kin$jday,
       pred_lower = apply(ppreds, 2, quantile, 0.025),
       pred_upper = apply(ppreds, 2, quantile, 0.975)) -> pred_obs_plot

# Make Figure S4
pred_obs_plot %>% 
  ggplot(aes(x = predicted, y = observed)) +
  geom_point(size = 3, alpha = 0.5) +
  geom_segment(aes(x = pred_lower, xend = pred_upper, y = observed, yend = observed), alpha = 0.2) +
  geom_abline(aes(intercept = 0, slope = 1), linewidth = 2, color = "dodgerblue") +
  ylim(95, 210) + xlim(95, 210) +
  annotate("text", label = bquote(R^2 == .(round(r2,3))), x = 100, y = 210, size = 7) +
  labs(x = "Predicted day of year", y = "Observed day of year") -> pred_obs

ppcs <- ppc_mean + ppc_sd

png("figs/FigS4_modelperform.png", height = 8.82, width = 10, res = 300, units = "in")
pred_obs / ppcs + plot_layout(heights = c(2,1)) + plot_annotation(tag_levels = "a",
                                                                  tag_prefix = "(",
                                                                  tag_suffix = ")")
dev.off()

## Create graphics - Results ####
# Preliminary results graphs
theme_set(theme_bw(base_size = 16))

# Get marginal means for PC x density x gravel interaction
pc1_plot <- sjPlot::plot_model(brms_lin_nokin, type = "emm", terms = c("pc1_sc", "density", "gravel"))
pc2_plot <- sjPlot::plot_model(brms_lin_nokin, type = "emm", terms = c("pc2_sc", "density", "gravel"))
# Get marginal means for site 
site_plot <- sjPlot::plot_model(brms_lin_nokin, type = "emm", terms = c("site"))

# Formatting and relabeling site level information
tibble(jday = site_plot$data$predicted,
       site = c("Cold\naseasonal (SS)", "Cool\nseasonal (CH)",
                "Hot\nseasonal (WI)", "Cool\naseasonal (BA)"),
       lower = site_plot$data$conf.low,
       upper = site_plot$data$conf.high) %>% 
  mutate(site = factor(site, levels = c("Cold\naseasonal (SS)",
                                        "Cool\nseasonal (CH)",
                                        "Cool\naseasonal (BA)",
                                        "Hot\nseasonal (WI)")))-> sum_stat_site

# Make Figure 3a
phen_flower_kin %>% 
  mutate(site = case_when(site == "SS" ~ "Cold\naseasonal (SS)",
                          site == "CH" ~ "Cool\nseasonal (CH)",
                          site == "WI" ~ "Hot\nseasonal (WI)",
                          site == "BA" ~ "Cool\naseasonal (BA)")) %>% 
  mutate(site = factor(site, levels = c("Cold\naseasonal (SS)",
                                        "Cool\nseasonal (CH)",
                                        "Cool\naseasonal (BA)",
                                        "Hot\nseasonal (WI)"))) %>% 
  ggplot(aes(x = site, y = jday)) +
  geom_jitter(height = 0, width = 0.3, aes(shape = site, color = site), size = 3, alpha = 0.1) +
  geom_segment(data = sum_stat_site, aes(x = site, y = lower, yend = upper, xend = site),
                  linewidth = 0.5) +
  geom_segment(data = sum_stat_site, aes(x = c(0.9, 1.9, 3.9, 2.9), y = jday,
                                         yend = jday, xend = c(1.1, 2.1, 4.1, 3.1)),
             linewidth = 0.5) +
  labs("") +
  scale_shape_manual(values = c(24,22,23,25)) +
  labs(y = "First day of flowering",
       x = "Site\n(cool & wet → hot & dry)") +
  theme(legend.position = "none") +
  scale_color_manual(values = c("#332288","#AA4499", "#44AA99", "#6699CC")) +
  scale_fill_manual(values = c("#332288","#AA4499", "#44AA99", "#6699CC")) -> site_subplot

# Create plotting dataset that is formatted 
phen_flower_kin_plot <- phen_flower_kin %>% 
  mutate(gravel = ifelse(gravel == "white", "High (black)", "Low (white)"),
         density = ifelse(density == "hi", "High", "Low"),
         site = case_when(site == "SS" ~ "Cold aseasonal (SS)",
                                 site == "CH" ~ "Cool seasonal (CH)",
                                 site == "WI" ~ "Hot seasonal (WI)",
                                 site == "BA" ~ "Cool aseasonal (BA)"))

# Create PC1 predictions part of Figure 3
tibble(pc = pc1_plot$data$x * sd_pc1 + mean_pc1,
       jday = pc1_plot$data$predicted,
       lower = pc1_plot$data$conf.low,
       upper = pc1_plot$data$conf.high,
       gravel = pc1_plot$data$facet,
       density = pc1_plot$data$group) %>% 
  mutate(gravel = ifelse(gravel == "white", "Low temp.\n(white gravel)", "High temp.\n(black gravel)"),
         density = ifelse(density == "hi", "High", "Low"),
         pc_type = "PC 1 (cool & wet → hot & dry)") -> pc1_dat

# Create PC2 predictions part of Figure 3
tibble(pc = pc2_plot$data$x * sd_pc2 + mean_pc2,
       jday = pc2_plot$data$predicted,
       lower = pc2_plot$data$conf.low,
       upper = pc2_plot$data$conf.high,
       gravel = pc2_plot$data$facet,
       density = pc2_plot$data$group) %>% 
  mutate(gravel = ifelse(gravel == "white", "Low temp.\n(white gravel)", "High temp.\n(black gravel)"),
         density = ifelse(density == "hi", "High", "Low"),
         pc_type = "PC 2 (low → high seasonality)") -> pc2_dat

# Create Figure 3b
rbind(pc1_dat, pc2_dat) %>% 
  ggplot(aes(pc, jday, color = density, fill = density)) +
  geom_line(aes(linetype = density), linewidth = 0.8) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = density), alpha = 0.4, color = NA) +
  facet_grid(gravel ~ pc_type, scales = "free_x") +
  labs(x = "Source climate PC axis value",
       y = "First day of flowering",
       color = "Density",
       fill = "Density",
       linetype = "Density") +
  scale_color_manual(values = c("#888888", "#DDCC77")) +
  scale_fill_manual(values = c("#888888", "#DDCC77")) +
  theme(legend.position = "top") -> pc_subplot

png("figs/Fig3_int.png", height = 9.5, width = 7.8, res = 300, units = "in")
site_subplot + pc_subplot + plot_annotation(tag_levels = "a", tag_prefix = "(",
                                            tag_suffix = ")") +
  plot_layout(heights = c(1,2)) 
dev.off()

# Create Fig S4 - density x gravel interaction across sites
pred_dat_int2 <- sjPlot::plot_model(brms_lin_nokin, terms = c("density", "gravel", "site"), type = "emm")
tibble(density = pred_dat_int2$data$x,
       gravel = pred_dat_int2$data$group,
       site = pred_dat_int2$data$facet,
       jday = pred_dat_int2$data$predicted,
       lower = pred_dat_int2$data$conf.low,
       upper = pred_dat_int2$data$conf.high) %>% 
  mutate(gravel = ifelse(gravel == "black", "High (black)", "Low (white)"),
         density = ifelse(density == 1, "High", "Low"),
         site = case_when(site == "SS" ~ "Cold aseasonal (SS)",
                          site == "CH" ~ "Cool seasonal (CH)",
                          site == "WI" ~ "Hot seasonal (WI)",
                          site == "BA" ~ "Cool aseasonal (BA)")) %>% 
  mutate(Site = factor(site, levels = c("Cold aseasonal (SS)",
                                        "Cool seasonal (CH)",
                                        "Cool aseasonal (BA)",
                                        "Hot seasonal (WI)"))) %>%
  ggplot(aes(x = gravel, y = jday, shape = site, color = density, fill = density)) +
  geom_jitter(data = phen_flower_kin_plot, aes(x = gravel, y = jday, shape = site), 
              position = position_jitterdodge(
                jitter.width = 0.2,
                jitter.height = 0,
                dodge.width = 0.75,
                seed = NA
              ), size = 0.8, alpha = 0.1) +
  geom_errorbar(aes(x = gravel, ymin = lower, ymax = upper),
                width = 0, position = position_dodge(width = 0.75), color = "black") +
  geom_point(aes(shape = site),size = 1.5, position = position_dodge(width = 0.75), color = "black") +
  facet_wrap(~site) +
  scale_color_manual(values = c("#888888", "#DDCC77")) +
  scale_shape_manual(values = c(24,23,22,25), guide = "none") +
  scale_fill_manual(values = c("#888888", "#DDCC77")) +
  labs(fill = "Density",
       color = "Density",
       x = "Temp. treatment (gravel)",
       y = "Day of year",
       shape = "") +
  guides(fill=guide_legend(override.aes=list(shape=21))) -> int_plot2
  
png("figs/FigS5_int.png", height = 9.5, width = 9.5, res = 300, units = "in") 
int_plot2
dev.off()

# Collect genotype-level information for genotype and GxE plots
as_draws_df(brms_lin_nokin) -> obj

# Collect random intercepts for each genotype
obj %>% 
  select(contains("r_genotype") & contains("Intercept")) %>% 
  gather(key = genotype, value = intercept) %>% 
  mutate(genotype = parse_number(genotype)) %>% 
  group_by(genotype) %>% 
  summarize(mean = mean(intercept)) %>% 
  arrange(genotype) %>% ungroup()-> random_intercepts

# Collect regression coefficients for PC1 and PC2 (genotype-level covariate)
phen_flower_kin %>% 
  group_by(genotype) %>% 
  summarize(pc1_sc = mean(pc1_sc),
            pc2_sc = mean(pc2_sc)) %>% 
  mutate(genotype = as.numeric(as.character(genotype))) %>% 
  arrange(genotype) %>% ungroup() -> pcs

# Get predicted mean of each genotype y = B0 + mu_genotype + B1*pc1_genotype + B2*pc2_genotype
pred <- mean(obj$b_Intercept) + random_intercepts$mean + pcs$pc1_sc * mean(obj$b_pc1_sc) + pcs$pc2_sc * mean(obj$b_pc2_sc)

# Join predictions with PC data and unscale PC values for plotting 
cbind(pcs, jday = pred) -> preds_by_genotype
preds_by_genotype$pc1 <- preds_by_genotype$pc1_sc * sd_pc1 + mean_pc1
preds_by_genotype$pc2 <- preds_by_genotype$pc2_sc * sd_pc2 + mean_pc2

# Create vector of PC1 for continuous plotting of response across PC1 axis
new_dat_pc1 <- seq(min(phen_flower_kin$pc1_sc), max(phen_flower_kin$pc1_sc), length.out = 10)
pred_pc1 <- mean(obj$b_Intercept) + mean(obj$b_pc1_sc)*new_dat_pc1 

# Create vector of PC2 for continuous plotting of response across PC2 axis
new_dat_pc2 <- seq(min(phen_flower_kin$pc2_sc), max(phen_flower_kin$pc2_sc), length.out = 10)
pred_pc2 <- mean(obj$b_Intercept) + mean(obj$b_pc2_sc)*new_dat_pc2

# Create storage and loop through predictions to get predicted means and
# quantiles for PC1
lower_pc1 <- NULL
upper_pc1 <- NULL
for (i in 1:10){
  temp <- obj$b_Intercept + obj$b_pc1_sc*new_dat_pc1[i]
  lower_pc1[i] <- quantile(temp, 0.025)
  upper_pc1[i] <- quantile(temp, 0.975)
}

# Create storage and loop through predictions to get predicted means and
# quantiles for PC2
lower_pc2 <- NULL
upper_pc2 <- NULL
for (i in 1:10){
  temp <- obj$b_Intercept + obj$b_pc2_sc*new_dat_pc2[i]
  lower_pc2[i] <- quantile(temp, 0.025)
  upper_pc2[i] <- quantile(temp, 0.975)
}

# Bring all PC1 prediction data together into one tibble (unscale + recenter)
tibble(pc1 = new_dat_pc1 * sd_pc1 + mean_pc1,
       jday = pred_pc1,
       lower = lower_pc1,
       upper = upper_pc1) -> sum_stat_pc1

# Bring all PC2 prediction data together into one tibble (unscale and recenter)
tibble(pc2 = new_dat_pc2 * sd_pc2 + mean_pc2,
       jday = pred_pc2,
       lower = lower_pc2,
       upper = upper_pc2) -> sum_stat_pc2

# Make PC1 plot (Fig 4a)
preds_by_genotype %>% 
  ggplot(aes(x = pc1, y = jday)) +
  geom_point(shape = 21, size = 3, alpha = 0.5, fill = "black", stroke = 1.5) +
  geom_line(data = sum_stat_pc1, color = "mediumpurple1", linewidth = 1.5) +
  geom_ribbon(data = sum_stat_pc1,
              aes(ymin = lower, ymax = upper, fill = density), alpha = 0.4, color = NA, fill = "mediumpurple1") +
  labs(x = "PC 1: cool & wet → hot & dry",
       y = "First day of flowering",
       fill = "") +
  ylim(139, 170) -> pc1_genotype_subplot

# Make PC2 plot (Fig 4b)
preds_by_genotype %>% 
  ggplot(aes(x = pc2, y = jday)) +
  geom_point(shape = 21, size = 3, alpha = 0.5, fill = "black", stroke = 1.5) +
  geom_line(data = sum_stat_pc2, color = "mediumpurple1", linewidth = 1.5) +
  geom_ribbon(data = sum_stat_pc2,
              aes(ymin = lower, ymax = upper, fill = density), alpha = 0.4, color = NA, fill = "mediumpurple1") +
  labs(x = "PC 2: low → high seasonality",
       y = "",
       fill = "") +
  ylim(139, 170) -> pc2_genotype_subplot

png("figs/Fig4_sourceclimate.png", height = 5, width = 9, units = "in", res = 300)
pc1_genotype_subplot + pc2_genotype_subplot +
  plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")") & theme(legend.position = "bottom") 
dev.off()

## Genotype by environment graph ####

# Pull out random slopes for gravel treatments
obj %>% 
  select(contains("r_genotype") & contains("gravel1")) %>% 
  gather(key = genotype, value = slope) %>% 
  mutate(genotype = parse_number(genotype)) %>% 
  group_by(genotype) %>% 
  summarize(effect_gravel = mean(slope)) %>% 
  arrange(genotype) %>% ungroup() -> random_slopes_gravel

# Predict black gravel mean for each genotype
mean(obj$b_Intercept) + mean(obj$b_gravel1) + 
      random_intercepts$mean + random_slopes_gravel$effect_gravel +
      mean(obj$b_pc1_sc)*pcs$pc1_sc + mean(obj$b_pc2_sc)*pcs$pc2_sc -> black_genotype

# Predict white gravel mean for each genotype
mean(obj$b_Intercept) - mean(obj$b_gravel1) +
  random_intercepts$mean - random_slopes_gravel$effect_gravel +
    mean(obj$b_pc1_sc)*pcs$pc1_sc + mean(obj$b_pc2_sc)*pcs$pc2_sc -> white_genotype

# Get Spearman's rank correlation for black and white gravel
spearman.ci(black_genotype, white_genotype)

# Make inset histogram of different slopes
tibble(gravel_diff = white_genotype - black_genotype) %>% 
  ggplot(aes(x = gravel_diff)) +
  geom_density(fill = "gray47") +
  geom_rug(length = unit(0.05, "npc")) +
  xlim(8.5, 12.2) +
  labs(x = "Gravel effect (days)\n[White - black gravel]",
       y = "Density") +
  theme_classic(base_size = 14) -> gravel_effect

# Pull out random slopes for density treatments
obj %>% 
  select(contains("r_genotype") & contains("density1")) %>% 
  gather(key = genotype, value = slope) %>% 
  mutate(genotype = parse_number(genotype)) %>% 
  group_by(genotype) %>% 
  summarize(effect_density = mean(slope)) %>% 
  arrange(genotype) %>% ungroup() -> random_slopes_density

# Predict high density for each genotype
mean(obj$b_Intercept) + mean(obj$b_density1) + 
      random_intercepts$mean + random_slopes_density$effect_density +
      mean(obj$b_pc1_sc)*pcs$pc1_sc + mean(obj$b_pc2_sc)*pcs$pc2_sc -> high_genotype
# Predict low density for each genotype
mean(obj$b_Intercept) - mean(obj$b_density1) +
      random_intercepts$mean - random_slopes_density$effect_density +
      mean(obj$b_pc1_sc)*pcs$pc1_sc + mean(obj$b_pc2_sc)*pcs$pc2_sc -> low_genotype

# Get Spearman's rank correlation for high and low gravel
spearman.ci(high_genotype, low_genotype)

# Make inset histogram of different slopes
tibble(density_diff = low_genotype - high_genotype) %>% 
  ggplot(aes(x = density_diff)) +
  geom_density(fill = "gray47") +
  geom_rug(length = unit(0.05, "npc")) +
  xlim(4.5,7) +
  labs(x = "Density effect (days)\n[Low - high density]",
       y = "Density") +
  theme_classic(base_size = 14) -> density_effect

# Pull random slopes for each site (can only do SS, BA, and WI based on how the
# model is parameterized -- will fit same model reordering factor to get CH)
obj %>% 
  select(contains("r_genotype")) %>% 
  select(contains("site1") | contains("site2") | contains("site3")) %>% 
  gather(key = genotype, value = slope) %>% 
  mutate(site = case_when(grepl("site1", genotype) ~ "SS",
                          grepl("site2", genotype) ~ "BA",
                          grepl("site3", genotype) ~ "WI")) %>% 
  mutate(genotype = parse_number(genotype)) %>% 
  group_by(genotype, site) %>% 
  summarize(effect_site = mean(slope)) %>% 
  spread(key = site, value = effect_site) %>% 
  arrange(genotype) %>% ungroup() -> random_slopes_site

# Predict SS for each genotype
mean(obj$b_Intercept) + mean(obj$b_site1) + 
  random_intercepts$mean + random_slopes_site$SS +
  mean(obj$b_pc1_sc)*pcs$pc1_sc + mean(obj$b_pc2_sc)*pcs$pc2_sc -> SS_genotype

# Predict BA for each genotype
mean(obj$b_Intercept) + mean(obj$b_site2) + 
  random_intercepts$mean + random_slopes_site$BA +
  mean(obj$b_pc1_sc)*pcs$pc1_sc + mean(obj$b_pc2_sc)*pcs$pc2_sc -> BA_genotype

# Predict WI for each genotype
mean(obj$b_Intercept) + mean(obj$b_site3) + 
  random_intercepts$mean + random_slopes_site$WI +
  mean(obj$b_pc1_sc)*pcs$pc1_sc + mean(obj$b_pc2_sc)*pcs$pc2_sc -> WI_genotype

# Relevel factor for second model to get CH value
phen_flower_kin$site <- factor(phen_flower_kin$site, levels = c("CH", "WI", "SS", "BA"))

# brms_lin_CHlevel <- brm(
#   jday ~ 0 + Intercept + density * gravel * pc1_sc + density * gravel * pc2_sc +
#     site * pc1_sc + site * pc2_sc + (1 + density + gravel + site || genotype) +
#     (1 | plot_unique),
#   data = phen_flower_kin,
#   data2 = list(Amat = kin),
#   family = gaussian(),
#   chains = 3, cores = 1, iter = 7500,
#   # Set seed for reproducibility
#   seed = 4685
# )

#write_rds(brms_lin_CHlevel, "~/Desktop/phenology_nokin_final_CH.rds")
brms_lin_CHlevel <- read_rds("outputs/phenology_nokin_final_CH.rds")

# Pull draws from model fit
obj_CH <- as_draws_df(brms_lin_CHlevel)

# Pull random intercepts in second model
obj_CH %>% 
  select(contains("r_genotype") & contains("Intercept")) %>% 
  gather(key = genotype, value = intercept) %>% 
  mutate(genotype = parse_number(genotype)) %>% 
  group_by(genotype) %>% 
  summarize(mean = mean(intercept)) %>% 
  arrange(genotype) %>% ungroup()-> random_intercepts_CH

# Pull random slopes for each site in second model
obj_CH %>% 
  select(contains("r_genotype") & contains("site1")) %>% 
  gather(key = genotype, value = slope) %>% 
  mutate(genotype = parse_number(genotype)) %>% 
  group_by(genotype) %>% 
  summarize(effect_CH = mean(slope)) %>% 
  arrange(genotype) %>% ungroup() -> random_slopes_CH

# Predict CH for each genotype
mean(obj_CH$b_Intercept) + mean(obj_CH$b_site1) + 
  random_intercepts_CH$mean + random_slopes_CH$effect_CH +
  mean(obj_CH$b_pc1_sc)*pcs$pc1_sc + mean(obj_CH$b_pc2_sc)*pcs$pc2_sc -> CH_genotype

# Get Spearman's rank correlation for all pairwise site combinations 
c(spearman.ci(SS_genotype, BA_genotype)$estimate,
spearman.ci(SS_genotype, WI_genotype)$estimate,
spearman.ci(SS_genotype, CH_genotype)$estimate,
spearman.ci(BA_genotype, WI_genotype)$estimate,
spearman.ci(BA_genotype, CH_genotype)$estimate,
spearman.ci(CH_genotype, WI_genotype)$estimate) -> site_rhos

mean(site_rhos)
sd(site_rhos) 

# Make inset histogram of different slopes
tibble(site_diff = SS_genotype - CH_genotype) %>% 
  ggplot(aes(x = site_diff)) +
  geom_density(fill = "gray47") +
  geom_rug(length = unit(0.05, "npc")) +
  xlim(-2.5,3) +
  labs(x = "Site effect (days)\n[Cold aseasonal - cool seasonal]",
       y = "Density") +
  theme_classic(base_size = 14) -> site_effect

# Create gravel GxE subplot
cbind(genotype = random_intercepts$genotype, white = white_genotype, black = black_genotype) %>% 
  as_tibble() %>% 
  gather(key = gravel, value = jday, white:black) %>% 
  mutate(gravel = ifelse(gravel == "black", "Black", "White")) %>% 
  ggplot(aes(x = gravel, y = jday, group = genotype)) + 
  geom_point(size = 3, pch = 21, alpha = 0.5, fill = "black") +
  geom_line(alpha = 0.5)+
  labs(y = "First day of flowering", x = "Gravel") -> gravel_gxe

# Create density GxE subplot
cbind(genotype = random_intercepts$genotype, high = high_genotype, low = low_genotype) %>% 
  as_tibble() %>% 
  gather(key = density, value = jday, high:low) %>% 
  mutate(density = ifelse(density == "high", "High", "Low")) %>% 
  ggplot(aes(x = density, y = jday, group = genotype)) + 
  geom_point(size = 3, pch = 21, alpha = 0.5, fill = "black") +
  geom_line(alpha = 0.5) +
  labs(y = "First day of flowering", x = "Density")-> density_gxe 

# Create site GxE subplot
cbind(genotype = random_intercepts$genotype,
      `Cold\naseasonal (SS)` = SS_genotype,
      `Cool\naseasonal (BA)` = BA_genotype,
      `Hot\nseasonal (WI)` = WI_genotype,
      `Cool\nseasonal (CH)` = CH_genotype) %>% 
  as_tibble() %>% 
  gather(key = site, value = jday, `Cold\naseasonal (SS)`:`Cool\nseasonal (CH)`) %>% 
  mutate(site = factor(site, levels = c("Cold\naseasonal (SS)",
                                        "Cool\nseasonal (CH)",
                                        "Cool\naseasonal (BA)",
                                        "Hot\nseasonal (WI)"))) %>% 
  ggplot(aes(x = site, y = jday, group = genotype)) + 
  geom_point(size = 3, pch = 21, fill = "black", alpha = 0.5) +
  geom_line(alpha = 0.5) +
  labs(y = "First day of flowering", x = "Site\n (cool & wet → hot & dry)") +
  scale_y_continuous(limits = c(120,175), breaks = seq(120,170,by=10))-> site_gxe

# Create design matrix for spacing out plots
design <- c(
  area(1, 1, 6, 3),
  area(1, 4, 3, 6),
  area(4, 4, 6, 6),
  area(7, 1, 7, 2),
  area(7, 3, 7, 4),
  area(7, 5, 7, 6)
)

# Visualize plot layout
# plot(design)

png("figs/Fig5_gxe.png", height = 8.73, width = 12.66, res = 300, units = "in")
(site_gxe | gravel_gxe / density_gxe) /
  (site_effect + gravel_effect + density_effect) + 
  plot_layout(heights = c(4,1)) +
  plot_annotation(tag_levels = "a",
                  tag_prefix = "(",
                  tag_suffix = ")") 
dev.off()

