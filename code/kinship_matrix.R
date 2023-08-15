# This code fits the linear model for flowering time

# Load libraries
library(tidyverse); library(mgcv); library(gratia); library(geomtextpath);
library(here); library(readr); library(brms); library(RcppCNPy); library(lme4);
library(patchwork)

# Source in compiled data for the model
source(here("supp_code", "compile_data.R"))

# Fit linear model without herbivory plants and with herbivory plants
phen_flower_kin %>% 
  filter(herbivory == "Y") %>% 
  group_by(site, density) %>% 
  summarize(n = n()) # Most of these are from Sheep Station

# Fit linear model with all data
mod_all <- lmer(jday ~ density*gravel*pc1 + density*gravel*pc2 + site*pc1 + site*pc2 +
                  (density+gravel+site|genotype) + (1|block_unique) + (1|plot_unique), data = phen_flower_kin)

car::Anova(mod_all)

# Drop all herbivory instances
phen_flower_kin %>% 
  filter(herbivory != "Y") -> phen_flower_kin_noherb

# Refit model
mod_noherb <- lmer(jday ~ density*gravel*pc1 + density*gravel*pc2 + site*pc1 + site*pc2 +
                     (density+gravel+site|genotype) + (1|block_unique) + (1|plot_unique), data = phen_flower_kin_noherb)

# Seems like the results are basically the same so fit the model with all of the
# data for now

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

gxe_density + gxe_gravel + gxe_site 

predicted_means %>% 
  group_by(genotype, pc1) %>% 
  summarize(jday_mean = mean(jday_pred)) %>% 
  ggplot(aes(x = reorder(genotype, jday_mean), y = jday_mean, fill = pc1)) +
  geom_point(aes(fill = pc1), shape = 21, size = 4) +
  scale_fill_distiller(palette = "PiYG") +
  coord_flip() +
  theme(axis.text.y = element_blank()) +
  ylab("jday") + xlab("genotype") -> genotype_pc1

predicted_means %>% 
  group_by(genotype, pc2) %>% 
  summarize(jday_mean = mean(jday_pred)) %>% 
  ggplot(aes(x = reorder(genotype, jday_mean), y = jday_mean, fill = pc2)) +
  geom_point(aes(fill = pc2), shape = 21, size = 4) +
  scale_fill_distiller(palette = "PuOr") +
  coord_flip() +
  theme(axis.text.y = element_blank()) +
  ylab("jday") + xlab("genotype") -> genotype_pc2

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


# Save this first run as a model object
# saveRDS(brms_m1, "supp_data/brms_flower.rds")
library(bayesplot)

ppreds <- posterior_predict(brms_m1, draws = 500)
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

# Calculate heritability
v_animal <- (VarCorr(brms_m1, summary = FALSE)$genotype$sd[,1])^2
v_r <- (VarCorr(brms_m1, summary = FALSE)$residual$sd)^2
h.bwt.1 <- as.mcmc(v_animal / (v_animal + v_r))
summary(h.bwt.1)
