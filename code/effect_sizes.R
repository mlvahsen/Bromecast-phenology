library(tidybayes)

as_draws_df(brms_m1) -> obj

# Store parameters to unscale later
mean_pc1 <- attr(scale(phen_flower_kin$pc1),"scaled:center")
sd_pc1 <- attr(scale(phen_flower_kin$pc1),"scaled:scale")
mean_pc2 <- attr(scale(phen_flower_kin$pc2),"scaled:center")
sd_pc2 <- attr(scale(phen_flower_kin$pc2),"scaled:scale")

# effect of gravel
quantile(obj$b_gravel1, c(0.5, 0.025, 0.975))
quantile(exp(obj$b_Intercept - obj$b_gravel1) - exp(obj$b_Intercept + obj$b_gravel1), c(0.5, 0.025, 0.975))
# effect of density
quantile(obj$b_density1, c(0.5, 0.025, 0.975))
quantile(exp(obj$b_Intercept - obj$b_density1) - exp(obj$b_Intercept + obj$b_density1), c(0.5, 0.025, 0.975))
# effect of density:gravel
quantile(obj$`b_density1:gravel1`, c(0.5, 0.025, 0.975))
# effect of site
quantile(obj$b_site3, c(0.5, 0.025, 0.975))
quantile(exp(obj$b_Intercept + obj$b_site2) - exp(obj$b_Intercept + obj$b_site3), c(0.5, 0.025, 0.975))
# effect of pc1
quantile(obj$b_pc1_sc, c(0.5, 0.025, 0.975))
# effect of pc2
quantile(obj$b_pc2_sc, c(0.5, 0.025, 0.975))
# effect of density:gravel:PC1
quantile(obj$`b_density1:gravel1:pc1_sc`, c(0.5, 0.025, 0.975))
# effect of genotype
quantile(obj$sd_genotype__Intercept, c(0.5, 0.025, 0.975))
# biggest effect of genotype
exp(obj$b_Intercept + obj$`r_genotype[34,Intercept]` + obj$b_pc1_sc*0.700233  + obj$b_pc2_sc*-0.1474011) -> earliest
exp(obj$b_Intercept + obj$`r_genotype[48,Intercept]` + obj$b_pc1_sc*-0.04574092 + obj$b_pc2_sc*1.048184) -> latest
quantile(latest - earliest, c(0.5, 0.025, 0.975))

quantile(exp(obj$b_Intercept + obj$b_site.BAMintercept) - exp(obj$b_Intercept + obj$b_site.WIMintercept), c(0.5, 0.025, 0.975))

# Supplemental plot showing interactions between PC 1 x site and PC 2 x site
site_pc1 <- sjPlot::plot_model(brms_m1, type = "emm", terms = c("pc1_sc", "site"))
site_pc2 <- sjPlot::plot_model(brms_m1, type = "emm", terms = c("pc2_sc", "site"))
  
tibble(pc1_sc = site_pc1$data$x,
       jday = site_pc1$data$predicted,
       site = site_pc1$data$group,
       lower = site_pc1$data$conf.low,
       upper = site_pc1$data$conf.high) %>% 
  mutate(pc1 = pc1_sc * sd_pc1 + mean_pc1) %>% 
  ggplot(aes(x = pc1, y = jday, color = site)) +
  geom_line() + 
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = site), alpha = 0.2, color = NA) +
  labs(x = "PC 1 (cool & wet → hot & dry)", y = "Day of year", color = "Site", fill = "Site") +
  ggtitle("") +
  scale_color_manual(values = c("#D55E00", "#009E73", "#CC79A7", "#0072B2")) +
  scale_fill_manual(values = c("#D55E00", "#009E73", "#CC79A7", "#0072B2")) -> site_pc1_plot

tibble(pc2_sc = site_pc2$data$x,
       jday = site_pc2$data$predicted,
       site = site_pc2$data$group,
       lower = site_pc2$data$conf.low,
       upper = site_pc2$data$conf.high) %>% 
  mutate(pc2 = pc2_sc * sd_pc2 + mean_pc2) %>% 
  ggplot(aes(x = pc2, y = jday, color = site)) +
  geom_line() + 
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = site), alpha = 0.2, color = NA) +
  labs(x = "PC 2 (high → low temperature seasonality)", y = "", color = "Site", fill = "Site") +
  ggtitle("") +
  scale_color_manual(values = c("#D55E00", "#009E73", "#CC79A7", "#0072B2")) +
  scale_fill_manual(values = c("#D55E00", "#009E73", "#CC79A7", "#0072B2"))-> site_pc2_plot

png("figs/FigS6_source_int.png", height = 5, width = 11, units = "in", res = 300)
site_pc1_plot + site_pc2_plot +
  plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")") +
  plot_layout(guides = "collect")
dev.off()