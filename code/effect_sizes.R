library(tidybayes)

as_draws_df(brms_lin) -> obj

# Store parameters to unscale later
mean_pc1 <- attr(scale(phen_flower_kin$pc1),"scaled:center")
sd_pc1 <- attr(scale(phen_flower_kin$pc1),"scaled:scale")
mean_pc2 <- attr(scale(phen_flower_kin$pc2),"scaled:center")
sd_pc2 <- attr(scale(phen_flower_kin$pc2),"scaled:scale")

# effect of gravel
quantile(obj$b_gravel1, c(0.5, 0.025, 0.975))
quantile((obj$b_Intercept - obj$b_gravel1) - (obj$b_Intercept + obj$b_gravel1), c(0.5, 0.025, 0.975))
# effect of density
quantile(obj$b_density1, c(0.5, 0.025, 0.975))
quantile((obj$b_Intercept - obj$b_density1) - (obj$b_Intercept + obj$b_density1), c(0.5, 0.025, 0.975))
# effect of density:gravel
quantile(obj$`b_density1:gravel1`, c(0.5, 0.025, 0.975))
# effect of site
quantile(obj$b_site3, c(0.5, 0.025, 0.975))
quantile((obj$b_Intercept + obj$b_site2) - (obj$b_Intercept + obj$b_site3), c(0.5, 0.025, 0.975))
# effect of pc1
quantile(obj$b_pc1_sc, c(0.5, 0.025, 0.975))
# effect of pc2
quantile(obj$b_pc2_sc, c(0.5, 0.025, 0.975))
# effect of density:gravel:PC1
quantile(obj$`b_density1:gravel1:pc1_sc`, c(0.5, 0.025, 0.975))
# effect of genotype
quantile(obj$sd_genotype__Intercept, c(0.5, 0.025, 0.975))
# biggest effect of genotype (earliest  = genotype 64; latest = genotype 6)
genotype64_pc1_sc <- unique(phen_flower_kin %>% filter(genotype == 64) %>% pull(pc1_sc))
genotype64_pc2_sc <- unique(phen_flower_kin %>% filter(genotype == 64) %>% pull(pc2_sc))
(obj$b_Intercept + obj$`r_genotype[64,Intercept]` + obj$b_pc1_sc*genotype64_pc1_sc  + obj$b_pc2_sc*genotype64_pc2_sc) -> earliest

genotype6_pc1_sc <- unique(phen_flower_kin %>% filter(genotype == 6) %>% pull(pc1_sc))
genotype6_pc2_sc <- unique(phen_flower_kin %>% filter(genotype == 6) %>% pull(pc2_sc))
(obj$b_Intercept + obj$`r_genotype[6,Intercept]` + obj$b_pc1_sc*genotype6_pc1_sc + obj$b_pc2_sc*genotype6_pc2_sc) -> latest
quantile(latest - earliest, c(0.5, 0.025, 0.975))

# Supplemental plot showing interactions between PC 1 x site and PC 2 x site
site_pc1 <- sjPlot::plot_model(brms_lin, type = "emm", terms = c("pc1_sc", "site"))
site_pc2 <- sjPlot::plot_model(brms_lin, type = "emm", terms = c("pc2_sc", "site"))
  
tibble(pc1_sc = site_pc1$data$x,
       jday = site_pc1$data$predicted,
       site = site_pc1$data$group,
       lower = site_pc1$data$conf.low,
       upper = site_pc1$data$conf.high) %>% 
  mutate(site = case_when(site == "SS" ~ "Cold aseasonal (SS)",
                   site == "CH" ~ "Cool seasonal (CH)",
                   site == "WI" ~ "Hot seasonal (WI)",
                   site == "BA" ~ "Cool aseasonal (BA)")) %>% 
  mutate(site = factor(site, levels = c("Cold aseasonal (SS)",
                                        "Cool seasonal (CH)",
                                        "Cool aseasonal (BA)",
                                        "Hot seasonal (WI)"))) %>%
  mutate(pc1 = pc1_sc * sd_pc1 + mean_pc1) %>% 
  ggplot(aes(x = pc1, y = jday, color = site)) +
  geom_line() + 
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = site), alpha = 0.2, color = NA) +
  labs(x = "PC 1 (cool & wet → hot & dry)", y = "Day of year", color = "Site", fill = "Site") +
  ggtitle("") +
  scale_color_manual(values = c("#009E73", "#0072B2", "#D55E00", "#CC79A7")) +
  scale_fill_manual(values = c("#009E73", "#0072B2", "#D55E00", "#CC79A7")) -> site_pc1_plot

tibble(pc2_sc = site_pc2$data$x,
       jday = site_pc2$data$predicted,
       site = site_pc2$data$group,
       lower = site_pc2$data$conf.low,
       upper = site_pc2$data$conf.high) %>% 
  mutate(site = case_when(site == "SS" ~ "Cold aseasonal (SS)",
                          site == "CH" ~ "Cool seasonal (CH)",
                          site == "WI" ~ "Hot seasonal (WI)",
                          site == "BA" ~ "Cool aseasonal (BA)")) %>% 
  mutate(site = factor(site, levels = c("Cold aseasonal (SS)",
                                        "Cool seasonal (CH)",
                                        "Cool aseasonal (BA)",
                                        "Hot seasonal (WI)"))) %>%
  mutate(pc2 = pc2_sc * sd_pc2 + mean_pc2) %>% 
  ggplot(aes(x = pc2, y = jday, color = site)) +
  geom_line() + 
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = site), alpha = 0.2, color = NA) +
  labs(x = "PC 2 (high → low temperature seasonality)", y = "", color = "Site", fill = "Site") +
  ggtitle("") +
  scale_color_manual(values = c("#009E73", "#0072B2", "#D55E00", "#CC79A7")) +
  scale_fill_manual(values = c("#009E73", "#0072B2", "#D55E00", "#CC79A7"))-> site_pc2_plot

png("figs/FigS6_source_int.png", height = 5, width = 11, units = "in", res = 300)
site_pc1_plot + site_pc2_plot +
  plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")") +
  plot_layout(guides = "collect")
dev.off()
