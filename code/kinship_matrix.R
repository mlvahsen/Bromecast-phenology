# This code fits the linear model for flowering time

# Load libraries
library(tidyverse); library(mgcv); library(gratia); library(geomtextpath);
library(here); library(readr); library(brms); library(RcppCNPy)

# Source in compiled data for the model
source(here("supp_code", "compile_data.R"))

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
sjPlot::plot_model(brms_m1, type = "emm", terms = c("pc2","density", "gravel")) +
  scale_color_manual(values = c("orange", "dodgerblue")) + ggtitle("")
dev.off()

png("figs/prelim_site.png", height = 5, width = 5.5, res = 300, units = "in")
sjPlot::plot_model(brms_m1, type = "emm", terms = c("gravel","density","site")) +
  scale_color_manual(values = c("orange", "dodgerblue")) + ggtitle("")
dev.off()
  

mod <- lmer(jday ~ density * gravel * pc1 + density * gravel * pc2 + site * pc1 + site * pc2 +
              (1|block_unique) + (1|plot_unique) + (1+density+gravel+site|genotype), data = phen_flower_kin)

anova(mod)
sjPlot::plot_model(mod, type = "pred", pred.type = "re", terms = c("site", "genotype"), se = F) +
  geom_line() + theme(legend.position = "none")

ranef(mod)


emmeans::emmeans(mod, ~site)

mod2 <- lmer(jday ~ density * gravel + site * pc1 + site * pc2 +
               (1|block_unique) + (1|plot_unique) + (1|genotype), data = phen_flower_kin)

anova(mod, mod2)

plot_model(mod, type = "pred", pred.type = "re", terms = c("density", "genotype"), se = F) +
  geom_line()

end <- Sys.time()

# Time difference of 1.263403 hours

# Save this first run as a model object
# saveRDS(brms_m1, "supp_data/brms_flower.rds")
brms_m1 <- readRDS("supp_data/brms_flower.rds")



newdata <- expand_grid(gravel = c("black", "white"),
                       density = c("hi", "lo"),
                       genotype = unique(phen_flower_kin$genotype),
                       site = unique(phen_flower_kin$site),
                       pc1 = mean(phen_flower_kin$pc1),
                       pc2 = mean(phen_flower_kin$pc2),
                       block_unique = "SS_1",
                       plot_unique = "SS_4_2")

predict(mod, newdata) -> preds

newdata %>% 
  group_by(density, genotype) %>% 
  summarize(mean = mean(preds)) %>% 
  ggplot(aes(x = density, y = mean, group = genotype)) +
  geom_line()

newdata$preds <- preds

epred_draws(brms_m1, newdata = newdata, re_formula = NULL, allow_new_levels = T) %>% 
  group_by(density, genotype) %>% 
  summarize(mean = mean(.epred)) %>% 
  spread(key = density, value = mean) %>% 
  mutate(diff = hi - lo) %>% 
  arrange(diff)
  

  

  
  ggplot(all_regions_autonomy_dist %>% filter(genotype %in% 1:10), 
         aes(x = .epred, y = genotype, 
             fill = site)) +
    stat_halfeye() +
    labs(x = "Predicted media index", y = NULL,
         fill = "Opposition parties allowed",
         subtitle = "Posterior predictions") +
    theme(legend.position = "bottom")

summary(brms_m1)

# Assess model fit
plot(brms_m1) # Looks pretty good!

# Calculate emmeans
emmeans::emmeans(brms_m1, ~density)

# Plot model
theme_set(theme_bw(base_size = 14))
sjPlot::plot_model(brms_m1, type = "emm", terms = c("density", "gravel", "site"))
sjPlot::plot_model(brms_m1, type = "emm", terms = c("pc1", "site"))
sjPlot::plot_model(brms_m1, type = "emm", terms = c("pc2", "site"))

ranef(brms_m1)

plot_model(brms_m1, type = "pred", 
           terms = c("density", "r_genotype [1:3]"), 
           pred.type = "re")

fitted <- fitted(brms_m1)
post_pred <- posterior_predict(brms_m1)
lower <- apply(post_pred, quantile, probs = c(0.025), MARGIN = 2)
upper <- apply(post_pred, quantile, probs = c(0.975), MARGIN = 2)

phen_flower_kin %>% 
  bind_cols(as_tibble(fitted))-> phen_fitted 
  

#group_by(genotype, density) %>% 
  #summarize(mean = mean(Estimate)) %>% 
phen_fitted %>% 
  ggplot(aes(x = jday, y = Estimate)) + geom_point() +
  geom_linerange(aes(ymin = lower, ymax = upper)) +
  geom_abline(aes(intercept = 0, slope = 1))

give_me_R2 <- function(preds,actual){
  rss <- sum(( preds - actual ) ^ 2)  ## residual sum of squares
  tss <- sum((actual - mean(actual)) ^ 2)  ## total sum of squares
  rsq <- 1 - rss/tss
  return(rsq)
}

random_persons <- sample(unique(phen_flower_kin$genotype), size = 20)
conditional_effects(brms_m1, type = "pred", effects = "density", 
                    re_formula = NULL, 
                    conditions = tibble(genotype = c("64", "48"))
)

brms_m1 %>% 
  emmeans(~density+genotype, re_formula = ~(1|genotype))

give_me_R2(phen_fitted$Estimate, phen_fitted$jday)

summary(brms_m1)

# Calculate heritability
v_animal <- (VarCorr(brms_m1, summary = FALSE)$genotype$sd[,1])^2
v_r <- (VarCorr(brms_m1, summary = FALSE)$residual$sd)^2
h.bwt.1 <- as.mcmc(v_animal / (v_animal + v_r))
summary(h.bwt.1)
