# Read in and format raw logger and snapshot temperature for common garden x
# gravel combinations

# Load libraries
library(tidyverse); library(lubridate)

# Read in snapshot temperatures
temp_m <- read_csv("data/micrositetemperature.csv")
# Read in logger data
temp <- read_csv("data/BCtemploggers.csv")

# Cut all data to be between Jan-1 and Jun-22 (where most sites x gravel have
# continuous data and encompasses most of growing season)
clean_dat <- temp %>% filter(Date >= "2022-01-01" & Date <= "2022-06-22")

# Set plotting theme
theme_set(theme_bw(base_size = 16))

# Clean and format logger data for plotting. Create Figure S3
clean_dat %>% 
  mutate(Site = case_when(Site == "Balzor" ~ "Cool aseasonal (BA)",
                          Site == "Cheyenne" ~ "Cool seasonal (CH)",
                          Site == "Sheep" ~ "Cold aseasonal (SS)",
                          Site == "Wildcat" ~ "Hot seasonal (WI)")) %>% 
  mutate(Site = factor(Site, levels = c("Cold aseasonal (SS)",
                                        "Cool seasonal (CH)",
                                        "Cool aseasonal (BA)",
                                        "Hot seasonal (WI)"))) %>% 
  mutate(Site = case_when(Site == "Balzor" ~ "Baltzor",
                          T ~ Site)) %>% 
  mutate(Color = case_when(Color == "Black" ~ "High (black)",
                           Color == "White" ~ "Low (white)")) %>% 
  ggplot(aes(x = Date, y = Temp_C, color = Color)) +
  geom_line(linewidth = 0.8) +
  facet_wrap(~Site) +
  scale_color_manual(values = c("black", "gray67")) +
  labs(y = "Temperature (°C) at 0-5 cm soil depth", x = "Date", color = "Temp. treatment\n(Gravel color)") -> soil05

png("figs/FigS3_SoilTemp5cm.png", width = 8, height = 6, res = 300, units = "in")
soil05
dev.off()

# Get mean and standard deviations for Table S1.
clean_dat %>% 
  select(Site, Date, Color, Temp_C) %>% 
  spread(key = Color, value = Temp_C) %>% 
  mutate(Diff = Black - White) %>% 
  group_by(Site) %>% 
  summarize(mean_black = mean(Black, na.rm = T),
            sd_black = sd(Black, na.rm = T)/sqrt(n()),
            mean_white = mean(White, na.rm = T),
            sd_white = sd(White, na.rm = T)/sqrt(n()),
            mean_diff = mean(Diff, na.rm = T),
            sd_diff = sd(Diff, na.rm = T)/sqrt(n()))

# Clean and format snapshot data for plotting. Create Figure S4
temp_m %>% 
  mutate(Site = case_when(Site == "Balzor" ~ "Cool aseasonal (BA)",
                          Site == "Cheyenne" ~ "Cool seasonal (CH)",
                          Site == "SheepStation" ~ "Cold aseasonal (SS)",
                          Site == "Wildcat" ~ "Hot seasonal (WI)")) %>% 
  mutate(Site = factor(Site, levels = c("Cold aseasonal (SS)",
                                        "Cool seasonal (CH)",
                                        "Cool aseasonal (BA)",
                                        "Hot seasonal (WI)"))) %>% 
  mutate(Color = case_when(Color == "Black" ~ "High (black)",
                           Color == "White" ~ "Low (white)")) %>% 
  mutate(Date = mdy(Date)) %>% 
  group_by(Site, Date, Color) %>% 
  summarize(mean = mean(SoilT_1cm, na.rm = T),
            sd = sd(SoilT_1cm)) %>% 
  ungroup() %>% 
  ggplot(aes(x = Date, y = mean)) +
  geom_point(size = 3, aes(shape = Site, fill = Color)) +
  scale_fill_manual(values = c("black", "gray67")) +
  facet_wrap(~Site) +
  labs(x = "Date", y = "Temperature (°C) at 0-1 cm soil depth",
       fill = "Temp. treatment\n(Gravel color)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_shape_manual(values = c(24,22,23,25), guide = "none") +
  guides(fill=guide_legend(override.aes=list(shape=21))) -> soil_temp

png("figs/FigS2_soil_temp.png", width = 8, height = 6, res = 300, units = "in")
soil_temp
dev.off()