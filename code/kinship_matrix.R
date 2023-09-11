# This code fits the linear model for flowering time

# Load libraries
library(tidyverse); library(mgcv); library(gratia); library(geomtextpath);
library(here); library(readr); library(brms); library(RcppCNPy); library(lme4);
library(patchwork); library(bayesplot); library(usmap)

# Source in compiled data for the model
source(here("supp_code", "compile_data.R"))

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

brms_m1 <- brm(
  jday ~ 1 + density * gravel * pc1 + density * gravel * pc2 +
    site * pc1 + site * pc2 + (1 + density + gravel + site || gr(genotype, cov = Amat)) +
    (1 | block_unique) + (1 | plot_unique),
  data = phen_flower_kin,
  data2 = list(Amat = kin),
  family = gaussian(),
  chains = 3, cores = 1, iter = 1000
)

end <- Sys.time()

end - start

summary(brms_m1)

## Create graphics - Model checking ####

brms_m1 <- read_rds("~/Downloads/brms_output.rds")

# Generate posterior predictive distribution
ppreds <- posterior_predict(brms_m1, draws = 500)
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
  ylim(105, 210) + xlim(105, 210) +
  annotate("text", label = bquote(R^2 == .(round(r2,3))), x = 117, y = 210, size = 7) -> pred_obs

png("figs/prelim_modelperform.png", height = 5.1, width = 5.75, res = 300, units = "in")
pred_obs
dev.off()


## Create graphics - Results ####
# Preliminary results graphs

theme_set(theme_bw(base_size = 16))

png("figs/prelim_int.png", height = 4, width = 8, res = 300, units = "in")
sjPlot::plot_model(brms_m1, type = "emm", terms = c("pc1","density", "gravel")) +
  scale_color_manual(values = c("orange", "dodgerblue")) + ggtitle("")
dev.off()

png("figs/prelim_site.png", height = 5, width = 5.5, res = 300, units = "in")
sjPlot::plot_model(brms_m1, type = "emm", terms = c("gravel","density","site")) +
  scale_color_manual(values = c("orange", "dodgerblue")) + ggtitle("")
dev.off()

# Create genotype graph

# Create new data data frame for predictions from model
expand_grid(density = c("lo", "hi"),
            gravel = c("black", "white"),
            plot_unique = unique(phen_flower_kin$plot_unique),
            genotype = factor(unique(phen_flower_kin$genotype))) %>% 
  mutate(block_unique = paste(str_split_i(plot_unique, "_", i = 1),
                              str_split_i(plot_unique, "_", i = 2), sep = "_"),
         site = str_split_i(plot_unique, "_", i = 1)) -> new_data

# Match up pc1 and pc2 values for each genotype from true data
phen_flower_kin %>% 
  ungroup() %>% 
  select(genotype, pc1, pc2) %>% 
  distinct() -> genotype_pcs

merge(new_data, genotype_pcs) -> new_data_pcs

predict(brms_m1, new_data_pcs) -> preds

# Bind predictions together with new data frame
cbind(new_data_pcs, jday_pred = preds[,1]) -> predicted_means

# Genotypes ranked by PC1 and PC2
predicted_means %>% 
  group_by(genotype, pc1) %>% 
  summarize(jday_mean = mean(jday_pred)) %>% 
  ggplot(aes(x = reorder(genotype, jday_mean), y = jday_mean, fill = pc1)) +
  geom_rug(length = unit(0.01, "npc"), alpha = 0.5, sides = "l") +
  geom_point(aes(fill = pc1), shape = 21, size = 4) +
  scale_fill_distiller(palette = "PiYG") +
  theme_classic(base_size = 14) +
  theme(axis.text.x = element_blank()) +
  labs(y = "julian day", x = "genotype", fill = "PC 1") +
  geom_curve(aes(x = 4, xend = 10, y = 143.5, yend = 141), curvature = 0.2, linewidth = 0.3) +
  geom_curve(aes(x = 37, xend = 43, y = 153.5, yend = 158), curvature = -0.2, linewidth = 0.3) +
  geom_text(aes(x = 16, y = 141), label = "Pahrump, NV", fontface = "italic", color = "maroon", size = 5) + 
  geom_text(aes(x = 51, y = 158), label = "Flathead Lake, MT", fontface = "italic", color = "darkgreen", size = 5) -> genotype_pc1 

predicted_means %>% 
  group_by(genotype, pc2) %>% 
  summarize(jday_mean = mean(jday_pred)) %>% 
  ggplot(aes(x = reorder(genotype, jday_mean), y = jday_mean, fill = pc2)) +
  geom_rug(length = unit(0.01, "npc"), alpha = 0.5, sides = "l") +
  geom_point(aes(fill = pc2), shape = 21, size = 4) +
  scale_fill_distiller(palette = "PuOr") +
  theme_classic(base_size = 14) +
  theme(axis.text.x = element_blank()) +
  ylab("julian day") + xlab("genotype") +
  labs(fill = "PC 2") +
  geom_curve(aes(x = 46, xend = 52, y = 154, yend = 160), curvature = -0.2, linewidth = 0.3) + 
  geom_curve(aes(x = 79, xend = 74, y = 159, yend = 152), curvature = -0.3, linewidth = 0.3) +
  geom_text(aes(x = 58, y = 160), label = "Badlands, SD", fontface = "italic", color = "brown", size = 5) +  
  geom_text(aes(x = 67.5, y = 152), label = "Tulameen, BC", fontface = "italic", color = "darkorchid4", size = 5) -> genotype_pc2

png("figs/Fig2_GenotypePC.png", height = 6, width = 11.8, units = "in", res = 300)
genotype_pc1 + genotype_pc2 + plot_layout(nrow = 2) + plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")")
dev.off()

predicted_means %>% 
  group_by(genotype, density) %>% 
  summarize(jday_mean = mean(jday_pred)) %>% 
  ggplot(aes(x = density, y = jday_mean, group = genotype)) +
  geom_line(alpha = 0.5) + geom_point(size = 3, alpha = 0.5) +
  labs(y = "jday")-> gxe_density

predicted_means %>% 
  group_by(genotype, gravel) %>% 
  summarize(jday_mean = mean(jday_pred)) %>% 
  ggplot(aes(x = gravel, y = jday_mean, group = genotype)) +
  geom_line(alpha = 0.5) + geom_point(size = 3, alpha = 0.5) +
  labs(y = "jday")-> gxe_gravel

predicted_means %>% 
  group_by(genotype, site) %>% 
  filter(site %in% c("CH", "SS", "BA")) %>% 
  summarize(jday_mean = mean(jday_pred)) %>% 
  ggplot(aes(x = site, y = jday_mean, group = genotype)) +
  geom_line(alpha = 0.5) + geom_point(size = 3, alpha = 0.5) +
  labs(y = "jday") -> gxe_site

predicted_means %>% 
  group_by(genotype, site) %>% 
  filter(site == "WI") %>% 
  summarize(jday_mean = mean(jday_pred)) %>% 
  arrange(jday_mean) %>%  
  ungroup() -> test_dat 
  

test_dat %>% 
  ggplot(aes(x = 1, y = jday_mean, color = jday_mean)) +
  geom_point(size = 2) +
  scale_color_distiller(palette = "Spectral") -> plot

layer_data(plot) %>% 
  mutate(genotype = test_dat$genotype) -> color_info

predicted_means %>% 
  group_by(genotype, site) %>% 
  filter(site == "BA") %>% 
  summarize(jday_mean = mean(jday_pred)) %>% 
  ungroup() %>% 
  merge(color_info %>% select(genotype, colour)) %>% 
  mutate(color_id = factor(1:92)) -> test_dat2 
  

test_dat2 %>% 
  ggplot(aes(x = 1, y = jday_mean, color = color_id)) +
  geom_point() +
  scale_color_manual(values = rainbow(92)) +
  theme(legend.position = "none")
  

gxe_density + gxe_gravel + gxe_site 



design <- c(
  area(1, 1, 6, 1),
  area(1, 2, 6, 2),
  area(1, 3, 2, 4),
  area(3, 3, 4, 4),
  area(5, 3, 6, 4)
)

plot(design)

png("figs/prelim_gxe.png", height = 9, width = 11, res = 300, units = "in")
genotype_pc1 + genotype_pc2 + gxe_density + gxe_gravel + gxe_site +
  plot_layout(design = design, guides = "collect") +
  plot_annotation(tag_levels = "a") & theme(legend.position = "left")
dev.off()

## Posthoc calculations of long- and short-term climate effects ####

# Get predicted mean jday for each cg site
predicted_means %>% 
  group_by(site) %>% 
  summarize(mean_jday = mean(jday_pred)) -> site_means

# Get predicted mean jday for each genotype
predicted_means %>% 
  group_by(genotype) %>% 
  summarize(mean_jday = mean(jday_pred)) %>% 
  filter(mean_jday == max(mean_jday) | mean_jday == min(mean_jday)) -> genotype_max_means

# Source in air temperature data
source("supp_code/prism_wrangle.R")

# Get temperature means for each common garden site
site_temp %>% 
  group_by(site_code) %>% 
  summarize(mean = mean(tmean)) -> site_temp_means

# Get temperature means for each collection site
collect_temp %>% 
  group_by(site_code) %>% 
  group_by(site_code) %>% 
  summarize(mean = mean(tmean)) %>% 
  filter(mean == max(mean) | mean == min(mean))-> collect_temp_means

# Bind together jday and temp data sets
site_ests <- merge(site_means %>% dplyr::select(site_code = site, mean_jday), site_temp_means)
collect_ests <- cbind(genotype_max_means, collect_temp_means) %>% dplyr::select(names(site_ests))

# Calculate relative phenological difference between CG sites and across
# collection sites
tibble(jday_diff = c(site_ests$mean_jday[1] - site_ests$mean_jday[4],
                     site_ests$mean_jday[2] - site_ests$mean_jday[4],
                     site_ests$mean_jday[3] - site_ests$mean_jday[4],
                     collect_ests$mean_jday[1] - collect_ests$mean_jday[2]),
       temp_diff = c(site_ests$mean[4] - site_ests$mean[1],
                     site_ests$mean[4] - site_ests$mean[2],
                     site_ests$mean[4] - site_ests$mean[3],
                     collect_ests$mean[2] - collect_ests$mean[1]),
       comp = c("short (BA vs WI)", "short (CH vs WI)", "short (SS vs WI)", "long (genotype source)")) %>% 
  mutate(rel_shift = jday_diff / temp_diff) -> rel_phen_diffs 
  
# Create relative phenological difference plot
rel_phen_diffs %>% 
  ggplot(aes(x =comp, y = rel_shift )) +
  geom_point(size = 8) + 
  geom_segment( aes(x=comp, xend=comp, y=0, yend=rel_shift), linewidth = 1.5) +
  xlab("") + ylab(expression(atop(paste("relative phenological shift "), paste("(", Delta," jday / ", Delta," mean temp)")))) +
  theme_bw(base_size = 18) +
  scale_x_discrete(labels = c("long \n (genotype source)", "short \n (BA vs WI)", "short \n (CH vs WI)", "short \n (SS vs WI)")) +
  coord_flip() -> rel_phen_plot

png("figs/prelim_relphenshift.png", height = 6.2, width = 7.3, res = 300, units = "in")
rel_phen_plot
dev.off()

## Calculations for in-text ####

# Calculate heritability -- still need to check how to do this when including
# random slopes
v_animal <- (VarCorr(brms_m1, summary = FALSE)$genotype$sd[,1])^2
v_r <- (VarCorr(brms_m1, summary = FALSE)$residual$sd)^2
h.bwt.1 <- as.mcmc(v_animal / (v_animal + v_r))
summary(h.bwt.1)
