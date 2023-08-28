predicted_means %>% 
  group_by(site) %>% 
  summarize(mean_jday = mean(jday_pred)) -> site_means

predicted_means %>% 
  group_by(genotype) %>% 
  summarize(mean_jday = mean(jday_pred)) %>% 
  filter(mean_jday == max(mean_jday) | mean_jday == min(mean_jday)) -> genotype_max_means

library(stringr)

tibble(jday_diff = c(160.4102-135.4704, 160.9054-135.4704, 160.9212-135.4704, 168.6375-140.2296),
       temp_diff = c(2.99, 0.96, 1.4, 12.9),
       comp = c("short (SS vs WI)", "short (CH vs WI)", "short (BA vs WI)", "long (genotype source)")) %>% 
  mutate(rel_shift = jday_diff / temp_diff) %>% 
  ggplot(aes(x =comp, y = rel_shift )) +
  geom_point(size = 8) + 
  geom_segment( aes(x=comp, xend=comp, y=0, yend=rel_shift), size = 1.5) +
  geom_text(label = c(8, 26, 18, 2), color = "white", fontface = 2) +
  xlab("") + ylab(expression(atop(paste("relative phenological shift "), paste("(", Delta," jday / ", Delta," mean temp)")))) +
  theme_bw(base_size = 18) +
  scale_x_discrete(labels = c("long \n (genotype source)", "short \n (BA vs WI)", "short \n (CH vs WI)", "short \n (SS vs WI)")) +
  coord_flip() -> rel_phen

png("figs/rel_phen_plot.png", height = 6.3, width = 7.4, res = 300, units = "in")
rel_phen
dev.off()

