# This code fits the linear model for flowering time

# Load libraries
library(tidyverse); library(mgcv); library(gratia); library(geomtextpath);
library(here); library(readr); library(brms); library(RcppCNPy); library(lme4);
library(patchwork); library(bayesplot); library(usmap); library(egg)

# Source in compiled data for the model
source(here("supp_code", "compile_data.R"))

# Center and scale continuous variables
phen_flower_kin$pc1_sc <- scale(phen_flower_kin$pc1)[,1]
phen_flower_kin$pc2_sc <- scale(phen_flower_kin$pc2)[,1]

# Store parameters to unscale later
mean_pc1 <- attr(scale(phen_flower_kin$pc1),"scaled:center")
sd_pc1 <- attr(scale(phen_flower_kin$pc1),"scaled:scale")
mean_pc2 <- attr(scale(phen_flower_kin$pc2),"scaled:center")
sd_pc2 <- attr(scale(phen_flower_kin$pc2),"scaled:scale")

# ## Check effect of herbivorized plants on analysis ####
# # Fit linear model without herbivory plants and with herbivory plants
# phen_flower_kin %>% 
#   filter(herbivory == "Y") %>% 
#   group_by(site, density) %>% 
#   summarize(n = n()) # Most of these are from Sheep Station
# 
# # Fit linear model with all data
# mod_all <- lmer(jday ~ density*gravel*pc1 + density*gravel*pc2 + site*pc1 + site*pc2 +
#                   (density+gravel+site|genotype) + (1|block_unique) + (1|plot_unique), data = phen_flower_kin)
# 
# emmeans::emmeans(mod_all, ~site)
# 
# # Drop all herbivory instances
# phen_flower_kin %>% 
#   filter(herbivory != "Y") -> phen_flower_kin_noherb
# 
# # Refit model
# mod_noherb <- lmer(jday ~ density*gravel*pc1 + density*gravel*pc2 + site*pc1 + site*pc2 +
#                      (density+gravel+site|genotype) + (1|block_unique) + (1|plot_unique), data = phen_flower_kin_noherb)
# 
# # Seems like the results are basically the same so fit the model with all of the
# # data for now
# 
## Fit Bayesian model ####
# Fit Bayesian linear model

start <- Sys.time()

# Fit with summed contrasts to more easily infer global means
options(contrasts = c("contr.sum", "contr.poly"))

brms_m1 <- brm(
  jday ~ 1 + density * gravel * pc1_sc + density * gravel * pc2_sc +
    site * pc1_sc + site * pc2_sc + (1 + density + gravel + site || gr(genotype, cov = Amat)) +
    (1 | block_unique) + (1 | plot_unique),
  data = phen_flower_kin,
  data2 = list(Amat = kin),
  family = poisson(link = "log"),
  chains = 3, cores = 1, iter = 15000,
  # Set seed for reproducibility
  seed = 1234
)

brms_NB <- brm(
  jday ~ 1 + density * gravel * pc1_sc + density * gravel * pc2_sc +
    site * pc1_sc + site * pc2_sc + (1 + density + gravel + site || gr(genotype, cov = Amat)) +
    (1 | block_unique) + (1 | plot_unique),
  data = phen_flower_kin,
  data2 = list(Amat = kin),
  family = negbinomial(link = "log", link_shape = "log"),
  chains = 1, cores = 1, iter = 1000,
  # Set seed for reproducibility
  seed = 1234
)

brms_lin <- brm(
  jday ~ 1 + density * gravel * pc1_sc + density * gravel * pc2_sc +
    site * pc1_sc + site * pc2_sc + (1 + density + gravel + site || gr(genotype, cov = Amat)) +
    (1 | plot_unique),
  data = phen_flower_kin,
  data2 = list(Amat = kin),
  family = gaussian(),
  chains = 1, cores = 1, iter = 5000,
  # Set seed for reproducibility
  seed = 4685
)

# write_rds(brms_m1, "~/Desktop/phenology_final_mod.rds")
brms_m1 <- read_rds("~/Desktop/phenology_final_mod.rds")

end <- Sys.time()

end - start
# Time difference of 0.6171344 hours

## Create graphics - Model checking ####

# Generate posterior predictive distribution
ppreds <- posterior_predict(brms_lin, draws = 100)
# Check against distribution of data, predicted mean, and predicted SD
ppc_dens_overlay(phen_flower_kin$jday, ppreds[1:50, ])
ppc_stat(phen_flower_kin$jday, ppreds, stat = "mean")
ppc_stat(phen_flower_kin$jday, ppreds, stat = "sd")

# Make a plot of predicted and observed for model
give_me_R2 <- function(preds,actual){
  rss <- sum(( preds - actual ) ^ 2)  ## residual sum of squares
  tss <- sum((actual - mean(actual)) ^ 2)  ## total sum of squares
  rsq <- 1 - rss/tss
  return(rsq)
}

r2 <- give_me_R2(colMeans(ppreds), phen_flower_kin$jday)

tibble(predicted = colMeans(ppreds),
       observed = phen_flower_kin$jday) %>% 
  ggplot(aes(x = predicted, y = observed)) +
  geom_point(size = 3, alpha = 0.5) +
  geom_abline(aes(intercept = 0, slope = 1), linewidth = 2, color = "dodgerblue") +
  ylim(108, 210) + xlim(108, 210) +
  annotate("text", label = bquote(R^2 == .(round(r2,3))), x = 115, y = 210, size = 7) +
  labs(x = "Predicted", y = "Observed") -> pred_obs

png("figs/FigS2_modelperform.png", height = 5.1, width = 5.75, res = 300, units = "in")
pred_obs
dev.off()

## Create graphics - Results ####
# Preliminary results graphs
brms_m1 <- brms_lin

theme_set(theme_bw(base_size = 16))

pc1_plot <- sjPlot::plot_model(brms_m1, type = "emm", terms = c("pc1_sc", "density", "gravel"))
pc2_plot <- sjPlot::plot_model(brms_m1, type = "emm", terms = c("pc2_sc", "density", "gravel"))
site_plot <- sjPlot::plot_model(brms_m1, type = "emm", terms = c("site"))

tibble(jday = site_plot$data$predicted,
       site = c("Sheep Station (SS)", "Cheyenne (CH)",
                "Wildcat (WI)", "Baltzor (BA)"),
       lower = site_plot$data$conf.low,
       upper = site_plot$data$conf.high) -> sum_stat_site

phen_flower_kin %>% 
  mutate(site = case_when(site == "SS" ~ "Sheep Station (SS)",
                          site == "CH" ~ "Cheyenne (CH)",
                          site == "WI" ~ "Wildcat (WI)",
                          site == "BA" ~ "Baltzor (BA)")) %>% 
  ggplot(aes(x = site, y = jday)) +
  geom_jitter(height = 0, width = 0.3, aes(shape = site, color = site), size = 3, alpha = 0.1) +
  geom_segment(data = sum_stat_site, aes(x = site, y = lower, yend = upper, xend = site),
                  linewidth = 1) +
  geom_point(data = sum_stat_site, aes(x = site, y = jday, shape = site, fill= site),
             size = 2.5) +
  labs("") +
  scale_shape_manual(values = 22:25) +
  labs(y = "Day of year",
       x = "Site") +
  theme(legend.position = "none", axis.text.x =element_text(angle = 20, hjust = 1)) +
  scale_color_manual(values = c("#009E73","#0072B2", "#D55E00", "#CC79A7")) +
  scale_fill_manual(values = c("#009E73","#0072B2", "#D55E00", "#CC79A7")) -> site_subplot

phen_flower_kin_plot <- phen_flower_kin %>% 
  mutate(gravel = ifelse(gravel == "white", "White gravel", "Black gravel"),
         density = ifelse(density == "hi", "High", "Low"),
         site = case_when(site == "SS" ~ "Sheep Station (SS)",
                                 site == "CH" ~ "Cheyenne (CH)",
                                 site == "WI" ~ "Wildcat (WI)",
                                 site == "BA" ~ "Baltzor (BA)"))

tibble(pc1 = pc1_plot$data$x * sd_pc1 + mean_pc1,
       jday = pc1_plot$data$predicted,
       lower = pc1_plot$data$conf.low,
       upper = pc1_plot$data$conf.high,
       gravel = pc1_plot$data$facet,
       density = pc1_plot$data$group) %>% 
  mutate(gravel = ifelse(gravel == "white", "White gravel", "Black gravel"),
         density = ifelse(density == "hi", "High", "Low")) %>% 
  ggplot(aes(x = pc1, y = jday, color = density)) +
  geom_point(data = phen_flower_kin_plot, aes(x = pc1, y = jday, color = density),
             shape = 1, alpha = 0.1) +
  geom_line(aes(linetype = density), linewidth = 1.2) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = density), alpha = 0.4, color = NA) +
  facet_wrap(~gravel) +
  labs(x = "PC 1 (cool & wet → hot & dry)",
       y = "Day of year",
       color = "Density",
       fill = "Density",
       linetype = "Density") +
  scale_color_manual(values = c("gray47", "maroon")) +
  scale_fill_manual(values = c("gray47", "maroon")) -> pc1_subplot
  
tibble(pc2 = pc2_plot$data$x * sd_pc2 + mean_pc2,
       jday = pc2_plot$data$predicted,
       lower = pc2_plot$data$conf.low,
       upper = pc2_plot$data$conf.high,
       gravel = pc2_plot$data$facet,
       density = pc2_plot$data$group) %>% 
  mutate(gravel = ifelse(gravel == "white", "White gravel", "Black gravel"),
         density = ifelse(density == "hi", "High", "Low")) %>% 
  ggplot(aes(x = pc2, y = jday, color = density)) +
  geom_point(data = phen_flower_kin_plot, aes(x = pc2, y = jday, color = density),
             shape = 1, alpha = 0.1) +
  geom_line(aes(linetype = density), linewidth = 1.2) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = density), alpha = 0.4, color = NA) +
  facet_wrap(~gravel) +
  labs(x = "PC 2 (low → high seasonality)",
       y = "Day of year",
       color = "Density",
       fill = "Density",
       linetype = "Density") +
  scale_color_manual(values = c("gray47", "maroon")) +
  scale_fill_manual(values = c("gray47", "maroon")) +
  theme(axis.title.x = element_text(margin = margin(t = -24, unit = "pt")))-> pc2_subplot


png("figs/Fig2_int.png", height = 7, width = 12, res = 300, units = "in")
site_subplot + (pc1_subplot / pc2_subplot) + plot_annotation(tag_levels = "a", tag_prefix = "(",
                                            tag_suffix = ")") +
  plot_layout(guides = "collect", widths = c(1,2))
dev.off()

# Create Fig S4 - density x gravel interaction across sites
pred_dat_int2 <- sjPlot::plot_model(brms_m1, terms = c("density", "gravel", "site"), type = "emm")
tibble(density = pred_dat_int2$data$x,
       gravel = pred_dat_int2$data$group,
       site = pred_dat_int2$data$facet,
       jday = pred_dat_int2$data$predicted,
       lower = pred_dat_int2$data$conf.low,
       upper = pred_dat_int2$data$conf.high) %>% 
  mutate(gravel = ifelse(gravel == "black", "Black gravel", "White gravel"),
         density = ifelse(density == 1, "High", "Low"),
         site = case_when(site == "SS" ~ "Sheep Station (SS)",
                          site == "CH" ~ "Cheyenne (CH)",
                          site == "WI" ~ "Wildcat (WI)",
                          site == "BA" ~ "Baltzor (BA)")) %>% 
  ggplot(aes(x = gravel, y = jday, color = density, fill = density)) +
  geom_jitter(data = phen_flower_kin_plot, aes(x = gravel, y = jday), shape = 1,
              position = position_jitterdodge(
                jitter.width = 0.2,
                jitter.height = 0,
                dodge.width = 0.75,
                seed = NA
              ), size = 0.8, alpha = 0.1) +
  geom_errorbar(aes(x = gravel, ymin = lower, ymax = upper),
                width = 0, position = position_dodge(width = 0.75), color = "black") +
  geom_point(size = 1.5, position = position_dodge(width = 0.75), pch = 21, color = "black") +
  facet_wrap(~site) +
  scale_color_manual(values = c("gray47", "maroon")) +
  scale_fill_manual(values = c("gray47", "maroon")) +
  labs(fill = "Density",
       color = "Density",
       x = "Gravel",
       y = "Day of year") -> int_plot2
  
png("figs/FigS4_int.png", height = 9.5, width = 9.5, res = 300, units = "in") 
int_plot2
dev.off()

# Create genotype graph
as_draws_df(brms_m1) -> obj

obj %>% 
  select(contains("r_genotype") & contains("Intercept")) %>% 
  gather(key = genotype, value = intercept) %>% 
  mutate(genotype = parse_number(genotype)) %>% 
  group_by(genotype) %>% 
  summarize(mean = mean(intercept)) %>% 
  arrange(genotype) %>% ungroup()-> random_intercepts

phen_flower_kin %>% 
  group_by(genotype) %>% 
  summarize(pc1_sc = mean(pc1_sc),
            pc2_sc = mean(pc2_sc)) %>% 
  mutate(genotype = as.numeric(as.character(genotype))) %>% 
  arrange(genotype) %>% ungroup() -> pcs

pred <- mean(obj$b_Intercept) + random_intercepts$mean + pcs$pc1_sc * mean(obj$b_pc1_sc) + pcs$pc2_sc * mean(obj$b_pc2_sc)

cbind(pcs, jday = pred) -> preds_by_genotype
preds_by_genotype$pc1 <- preds_by_genotype$pc1_sc * sd_pc1 + mean_pc1
preds_by_genotype$pc2 <- preds_by_genotype$pc2_sc * sd_pc2 + mean_pc2

new_dat_pc1 <- seq(min(phen_flower_kin$pc1_sc), max(phen_flower_kin$pc1_sc), length.out = 10)
pred_pc1 <- mean(obj$b_Intercept) + mean(obj$b_pc1_sc)*new_dat_pc1 

new_dat_pc2 <- seq(min(phen_flower_kin$pc2_sc), max(phen_flower_kin$pc2_sc), length.out = 10)
pred_pc2 <- mean(obj$b_Intercept) + mean(obj$b_pc2_sc)*new_dat_pc2

lower_pc1 <- NULL
upper_pc1 <- NULL
for (i in 1:10){
  temp <- obj$b_Intercept + obj$b_pc1_sc*new_dat_pc1[i]
  lower_pc1[i] <- quantile(temp, 0.025)
  upper_pc1[i] <- quantile(temp, 0.975)
}

lower_pc2 <- NULL
upper_pc2 <- NULL
for (i in 1:10){
  temp <- obj$b_Intercept + obj$b_pc2_sc*new_dat_pc2[i]
  lower_pc2[i] <- quantile(temp, 0.025)
  upper_pc2[i] <- quantile(temp, 0.975)
}

tibble(pc1 = new_dat_pc1 * sd_pc1 + mean_pc1,
       jday = pred_pc1,
       lower = lower_pc1,
       upper = upper_pc1) -> sum_stat_pc1

tibble(pc2 = new_dat_pc2 * sd_pc2 + mean_pc2,
       jday = pred_pc2,
       lower = lower_pc2,
       upper = upper_pc2) -> sum_stat_pc2

preds_by_genotype %>% 
  ggplot(aes(x = pc1, y = jday)) +
  geom_point(shape = 21, size = 3, alpha = 0.7, fill = "black", stroke = 1.5) +
  geom_line(data = sum_stat_pc1, color = "black") +
  geom_line(data = sum_stat_pc1, aes(y = lower), color = "black", linetype = "dashed") +
  geom_line(data = sum_stat_pc1, aes(y = upper), color = "black", linetype = "dashed") +
  labs(x = "PC 1: cool & wet → hot & dry",
       y = "Day of year",
       fill = "") +
  ylim(140, 175) -> pc1_genotype_subplot

preds_by_genotype %>% 
  ggplot(aes(x = pc2, y = jday)) +
  geom_point(shape = 21, size = 3, alpha = 0.7, fill = "black", stroke = 1.5) +
  geom_line(data = sum_stat_pc2) +
  geom_line(data = sum_stat_pc2, aes(y = lower), color = "black", linetype = "dashed") +
  geom_line(data = sum_stat_pc2, aes(y = upper), color = "black", linetype = "dashed") +
  labs(x = "PC 2: low → high seasonality",
       y = "",
       fill = "") +
  ylim(140, 175) -> pc2_genotype_subplot

png("figs/Fig3_sourceclimate.png", height = 5, width = 9, units = "in", res = 300)
pc1_genotype_subplot + pc2_genotype_subplot +
  plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")") & theme(legend.position = "bottom") 
dev.off()

## Genotype by environment graph ####
obj %>% 
  select(contains("r_genotype") & contains("Intercept")) %>% 
  gather(key = genotype, value = intercept) %>% 
  mutate(genotype = parse_number(genotype)) %>% 
  group_by(genotype) %>% 
  summarize(mean = mean(intercept)) %>% 
  arrange(genotype) %>% ungroup()-> random_intercepts

obj %>% 
  select(contains("r_genotype") & contains("gravel1")) %>% 
  gather(key = genotype, value = slope) %>% 
  mutate(genotype = parse_number(genotype)) %>% 
  group_by(genotype) %>% 
  summarize(effect_gravel = mean(slope)) %>% 
  arrange(genotype) %>% ungroup() -> random_slopes_gravel

mean(obj$b_Intercept) + mean(obj$b_gravel1) + 
      random_intercepts$mean + random_slopes_gravel$effect_gravel +
      mean(obj$b_pc1_sc)*pcs$pc1_sc + mean(obj$b_pc2_sc)*pcs$pc2_sc -> black_genotype
mean(obj$b_Intercept) - mean(obj$b_gravel1) +
  random_intercepts$mean - random_slopes_gravel$effect_gravel +
    mean(obj$b_pc1_sc)*pcs$pc1_sc + mean(obj$b_pc2_sc)*pcs$pc2_sc -> white_genotype

obj %>% 
  select(contains("r_genotype") & contains("density1")) %>% 
  gather(key = genotype, value = slope) %>% 
  mutate(genotype = parse_number(genotype)) %>% 
  group_by(genotype) %>% 
  summarize(effect_density = mean(slope)) %>% 
  arrange(genotype) %>% ungroup() -> random_slopes_density

mean(obj$b_Intercept) + mean(obj$b_density1) + 
      random_intercepts$mean + random_slopes_density$effect_density +
      mean(obj$b_pc1_sc)*pcs$pc1_sc + mean(obj$b_pc2_sc)*pcs$pc2_sc -> high_genotype
mean(obj$b_Intercept) - mean(obj$b_density1) +
      random_intercepts$mean - random_slopes_density$effect_density +
      mean(obj$b_pc1_sc)*pcs$pc1_sc + mean(obj$b_pc2_sc)*pcs$pc2_sc -> low_genotype

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

mean(obj$b_Intercept) + mean(obj$b_site1) + 
  random_intercepts$mean + random_slopes_site$SS +
  mean(obj$b_pc1_sc)*pcs$pc1_sc + mean(obj$b_pc2_sc)*pcs$pc2_sc -> SS_genotype

mean(obj$b_Intercept) + mean(obj$b_site2) + 
  random_intercepts$mean + random_slopes_site$BA +
  mean(obj$b_pc1_sc)*pcs$pc1_sc + mean(obj$b_pc2_sc)*pcs$pc2_sc -> BA_genotype

mean(obj$b_Intercept) + mean(obj$b_site3) + 
  random_intercepts$mean + random_slopes_site$WI +
  mean(obj$b_pc1_sc)*pcs$pc1_sc + mean(obj$b_pc2_sc)*pcs$pc2_sc -> WI_genotype

phen_flower_kin$site <- factor(phen_flower_kin$site, levels = c("CH", "WI", "SS", "BA"))

brms_lin_CHlevel <- brm(
  jday ~ 1 + density * gravel * pc1_sc + density * gravel * pc2_sc +
    site * pc1_sc + site * pc2_sc + (1 + density + gravel + site || gr(genotype, cov = Amat)) +
    (1 | plot_unique),
  data = phen_flower_kin,
  data2 = list(Amat = kin),
  family = gaussian(),
  chains = 1, cores = 1, iter = 5000,
  # Set seed for reproducibility
  seed = 4685
)

obj_CH <- as_draws_df(brms_lin_CHlevel)

obj_CH %>% 
  select(contains("r_genotype") & contains("Intercept")) %>% 
  gather(key = genotype, value = intercept) %>% 
  mutate(genotype = parse_number(genotype)) %>% 
  group_by(genotype) %>% 
  summarize(mean = mean(intercept)) %>% 
  arrange(genotype) %>% ungroup()-> random_intercepts_CH

obj_CH %>% 
  select(contains("r_genotype") & contains("site1")) %>% 
  gather(key = genotype, value = slope) %>% 
  mutate(genotype = parse_number(genotype)) %>% 
  group_by(genotype) %>% 
  summarize(effect_CH = mean(slope)) %>% 
  arrange(genotype) %>% ungroup() -> random_slopes_CH

mean(obj_CH$b_Intercept) + mean(obj_CH$b_site1) + 
  random_intercepts_CH$mean + random_slopes_CH$effect_CH +
  mean(obj_CH$b_pc1_sc)*pcs$pc1_sc + mean(obj_CH$b_pc2_sc)*pcs$pc2_sc -> CH_genotype

# Gravel gxe plot
cbind(genotype = random_intercepts$genotype, white = white_genotype, black = black_genotype) %>% 
  as_tibble() %>% 
  gather(key = gravel, value = jday, white:black) %>% 
  mutate(gravel = ifelse(gravel == "black", "Black", "White")) %>% 
  ggplot(aes(x = gravel, y = jday, group = genotype)) + 
  geom_point(size = 3, pch = 21, alpha = 0.5, fill = "black") +
  geom_line(alpha = 0.5) +
  ylim(135, 175) +
  labs(y = "Day of year", x = "Gravel") -> gravel_gxe

# Density gxe plot
cbind(genotype = random_intercepts$genotype, high = high_genotype, low = low_genotype) %>% 
  as_tibble() %>% 
  gather(key = density, value = jday, high:low) %>% 
  mutate(density = ifelse(density == "high", "High", "Low")) %>% 
  ggplot(aes(x = density, y = jday, group = genotype)) + 
  geom_point(size = 3, pch = 21, alpha = 0.5, fill = "black") +
  geom_line(alpha = 0.5) +
  ylim(135, 175) +
  labs(y = "Day of year", x = "Density")-> density_gxe 

# Site gxe plot
cbind(genotype = random_intercepts$genotype,
      `Sheep Station (SS)` = SS_genotype,
      `Baltzor (BA)` = BA_genotype,
      `Wildcat (WI)` = WI_genotype,
      `Cheyenne (CH)` = CH_genotype) %>% 
  as_tibble() %>% 
  gather(key = site, value = jday, `Sheep Station (SS)`:`Cheyenne (CH)`) %>% 
  ggplot(aes(x = site, y = jday, group = genotype)) + 
  geom_point(size = 3, pch = 21, fill = "black", alpha = 0.5) +
  geom_line(alpha = 0.5) +
  labs(y = "Day of year", x = "Site") + 
  ylim(120,175) +
  theme(axis.text.x = element_text(angle = 10, hjust = 1)) -> site_gxe

design <- c(
  area(1, 1, 4, 3),
  area(1, 4, 2, 5),
  area(3, 4, 4, 5)
)


png("figs/Fig4_gxe.png", height = 7, width = 12.4, res = 300, units = "in")
site_gxe + gravel_gxe + density_gxe +
  plot_layout(design = design, guides = "collect") +
  plot_annotation(tag_levels = "a",
                  tag_prefix = "(",
                  tag_suffix = ")") 
dev.off()
