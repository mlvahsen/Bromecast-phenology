library(tidyverse)

temp_m <- read_csv("data/micrositetemperature.csv")
temp <- read_csv("data/BCtemploggers.csv")

theme_set(theme_bw(base_size = 16))

temp %>%
  filter(Site == "Cheyenne" & Date < "2022-05-06") -> Cheyenne
temp %>% 
  filter(Site == "Wildcat" & Date > "2021-09-30" & Date < "2022-07-24") -> Wildcat
temp %>% 
  filter(Site == "Balzor" & Date < "2022-07-04") -> Baltzor
temp %>% 
  filter(Site == "Sheep" & Date < "2022-07-24" & Date > "2021-11-18") -> Sheep

clean_dat <- rbind(Cheyenne, Wildcat, Baltzor, Sheep)

clean_dat %>% 
  mutate(Site = case_when(Site == "Balzor" ~ "Baltzor",
                          T ~ Site)) %>% 
  ggplot(aes(x = Date, y = Temp_C, color = Color)) +
  geom_line(linewidth = 0.8) +
  facet_wrap(~Site) +
  scale_color_manual(values = c("black", "gray67")) +
  labs(y = "Temperature (°C) at 0-5 cm soil depth", x = "Date", color = "Gravel color") -> soil05

png("figs/FigS3.png", width = 8, height = 6, res = 300, units = "in")
soil05
dev.off()

clean_dat %>% 
  spread(key = Color, value = Temp_C) %>% 
  mutate(Diff = Black - White) %>% 
  group_by(Site) %>% 
  summarize(mean_black = mean(Black, na.rm = T),
            sd_black = sd(Black, na.rm = T)/sqrt(n()),
            mean_white = mean(White, na.rm = T),
            sd_white = sd(White, na.rm = T)/sqrt(n()),
            mean_diff = mean(Diff, na.rm = T))

clean_dat %>% 
  group_by(Site, Color) %>% 
  summarize(n = n())



temp_m %>% 
  ggplot(aes(x = Date, y = SoilT_1cm, color = Color, linetype = Density)) +
  geom_line(linewidth = 0.8) +
  facet_wrap(~Site) +
  scale_color_manual(values = c("black", "gray67")) +
  labs(y = "temperature (°C)", x = "date", color = "gravel color") +
  ggtitle("Soil temperature at 1 cm depth (°C)")

temp_m %>% 
  mutate(Site = case_when(Site == "Balzor" ~ "Baltzor",
                          T ~ Site))  %>% 
  group_by(Site, Date, Color) %>% 
  summarize(mean = mean(SoilT_1cm, na.rm = T),
            se = sd(SoilT_1cm)/sqrt(n())) %>% 
  ungroup() %>% 
  ggplot(aes(x = Date, y = mean, color = Color)) +
  geom_point(size = 3) +
  facet_wrap(~Site) +
  geom_errorbar(aes(x = Date, ymin = mean - se, ymax = mean + se), width = 0.1) +
  labs(x = "Date", y = "Temperature (°C) at 0-1 cm soil depth",
       color = "Gravel color") +
  scale_color_manual(values = c("black", "gray67")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) -> soil_temp

temp_m %>% 
  group_by(Site, Color) %>% 
  summarize(mean = mean(SoilT_1cm, na.rm = T)) %>% 
  spread(key = Color, value = mean) %>% 
  mutate(diff = Black - White)


png("figs/FigS2_soil_temp.png", width = 8, height = 6, res = 300, units = "in")
soil_temp
dev.off()
