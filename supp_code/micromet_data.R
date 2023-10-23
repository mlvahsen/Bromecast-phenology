library(tidyverse)

temp_m <- read_csv("data/micrositetemperature.csv")
temp <- read_csv("data/BCtemploggers.csv")

theme_set(theme_bw(base_size = 16))

temp %>% 
  ggplot(aes(x = Date, y = Temp_C, color = Color)) +
  geom_line(linewidth = 0.8) +
  facet_wrap(~Site) +
  scale_color_manual(values = c("black", "gray67")) +
  labs(y = "temperature (°C)", x = "date", color = "gravel color") +
  ggtitle("Microsite temperature (°C)") -> air_temp

png("figs/air_temp.png", width = 14.75, height = 9, res = 300, units = "in")
air_temp
dev.off()

temp_m %>% 
  ggplot(aes(x = Date, y = SoilT_1cm, color = Color, linetype = Density)) +
  geom_line(linewidth = 0.8) +
  facet_wrap(~Site) +
  scale_color_manual(values = c("black", "gray67")) +
  labs(y = "temperature (°C)", x = "date", color = "gravel color") +
  ggtitle("Soil temperature at 1 cm depth (°C)")

temp_m %>% 
  group_by(Site, Date, Color, Density) %>% 
  summarize(mean = mean(SoilT_1cm, na.rm = T),
            se = sd(SoilT_1cm)) %>% 
  ungroup() %>% 
  ggplot(aes(x = Date, y = mean, color = Color, shape = Density)) +
  geom_point(position = position_dodge(width = 0.5), size = 3) +
  facet_wrap(~Site) +
  geom_errorbar(aes(x = Date, ymin = mean - se, ymax = mean + se), width = 0.1,
                position = position_dodge(width = 0.5)) +
  labs(x = "date", y = "temperature (°C)", shape = "density treatment",
       color = "gravel color") +
  scale_color_manual(values = c("black", "gray67")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Soil temperature at 1 cm depth (°C)") -> soil_temp

png("figs/soil_temp.png", width = 14.75, height = 9, res = 300, units = "in")
soil_temp
dev.off()
