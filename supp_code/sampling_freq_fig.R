# Create a plot to depict when plants where surveyed for each site

## Load libraries ####
library(tidyverse); library(mgcv); library(gratia); library(geomtextpath);
library(here); library(readr); library(brms); library(RcppCNPy)

## Source code for genotype climate of origin ####
source(here("supp_code/climate_of_origin.R"))

## Read in Sheep Station data ####

# Read in derived phenology data
phen_SS <- read_csv("~/Git/Bromecast/gardens/deriveddata/SS2022_growthphenology_with_harvest.csv")
# Read in plant ID info
ids_SS <- read_csv("~/Git/Bromecast/gardens/deriveddata/SS2022_plantID.csv")
# Read in flagging data
flags_SS <- read_csv("~/Git/Bromecast/gardens/deriveddata/SS2022_flags.csv")
gardens <- read_csv("~/Git/Bromecast/gardens/rawdata/garden_treatments.csv")
# Merge together datasets
phen_id_SS <- merge(phen_SS, ids_SS)

# Rename 'garden' column and remove cum_plot column
gardens %>% 
  mutate(site = garden) %>% 
  dplyr::select(-cum_plot, -garden) -> gardens_sub

# Merge together datasets
phen_id_garden_SS <- merge(phen_id_SS, gardens_sub)

# Merge together datasets
phen_SS <- merge(phen_id_garden_SS, flags_SS)
phen_SS <- merge(phen_SS, genotype_PCclimate)

# Set appropriate factors for variables
phen_SS %>% 
  mutate(block = as.factor(block),
         plot = as.factor(plot),
         growout = as.factor(growout),
         density = as.factor(density),
         gravel = as.factor(gravel),
         site = as.factor(site),
         genotype = as.factor(genotype)) %>% 
  mutate(plot_unique = as.factor(paste(site, block, plot, sep = "_")),
         block_unique = as.factor(paste(site, block, sep = "_")))-> phen_SS


## Read in Boise data ####
# Read in derived phenology data
phen_Boise <- read_csv("~/Git/Bromecast/gardens/deriveddata/Boise2022_growthphenology_by_plantID.csv")
# Read in plant ID info
ids_Boise <- read_csv("~/Git/Bromecast/gardens/deriveddata/Boise2022_plantID.csv")
# Read in flagging data
flags_Boise <- read_csv("~/Git/Bromecast/gardens/deriveddata/Boise2022_flags.csv")

# Merge together datasets
phen_id_Boise <- merge(phen_Boise, ids_Boise)

# Merge together datasets
phen_id_garden_Boise <- merge(phen_id_Boise, gardens_sub)

# Merge together datasets
phen_Boise <- merge(phen_id_garden_Boise, flags_Boise)
phen_Boise <- merge(phen_Boise, genotype_PCclimate)

# Set appropriate factors for variables
phen_Boise %>% 
  mutate(block = as.factor(block),
         plot = as.factor(plot),
         growout = NA,
         density = as.factor(density),
         gravel = as.factor(gravel),
         site = as.factor(site),
         genotype = as.factor(genotype)) %>% 
  mutate(plot_unique = as.factor(paste(site, block, plot, sep = "_")),
         block_unique = as.factor(paste(site, block, sep = "_")))-> phen_Boise


## Read in Cheyenne data ####
# Read in derived phenology data
phen_CH <- read_csv("~/Git/Bromecast/gardens/deriveddata/CH2022_growthphenology_by_plantID.csv")
# Read in plant ID info
ids_CH <- read_csv("~/Git/Bromecast/gardens/deriveddata/CH2022_plantID.csv")
# Read in flagging data
flags_CH <- read_csv("~/Git/Bromecast/gardens/deriveddata/CH2022_flags.csv")

# Merge together datasets
phen_id_CH <- merge(phen_CH, ids_CH)

# Merge together datasets
phen_CH <- merge(phen_id_CH, flags_CH)
phen_CH <- merge(phen_CH, genotype_PCclimate)

# Set appropriate factors for variables
phen_CH %>% 
  mutate(block = as.factor(block),
         plot = as.factor(plot),
         growout = NA,
         density = as.factor(case_when(density == "high" ~ "hi",
                                       density == "low" ~ "lo")),
         gravel = as.factor(gravel),
         site = as.factor(site),
         genotype = as.factor(genotype)) %>% 
  mutate(plot_unique = as.factor(paste(site, block, plot, sep = "_")),
         block_unique = as.factor(paste(site, block, sep = "_")))-> phen_CH

## Merge datasets together ####
phen <- rbind(phen_SS %>% dplyr::select(-tillers), phen_Boise, phen_CH)

phen %>%
  filter(v %in% c("FG", "FP", "FB", "FX")) %>%
  group_by(plantID) %>%
  # Gets minimum day of flowering
  slice(which.min(jday)) -> phen_flower

phen_flower %>% 
  group_by(site) %>% 
  slice(which.min(jday)) %>% 
  select(site, jday) -> first_flower_site

phen %>% filter(site == "SS" & jday > 50) %>% pull(jday) %>% range(na.rm = T)
phen %>% filter(site == "WI" & jday > 50) %>% pull(jday) %>% range(na.rm = T)
phen %>% filter(site == "BA" & jday > 50) %>% pull(jday) %>% range(na.rm = T)
phen %>% filter(site == "CH" & jday > 50) %>% pull(jday) %>% range(na.rm = T)

colors <- c("#44AA99", "#332288", "#6699CC", "#AA4499")

phen %>% 
  group_by(jday, site) %>% 
  summarize(n = n()) %>% 
  filter(jday > 50) %>% 
  mutate(y = case_when(site == "SS" ~ 4,
                       site == "CH" ~ 3,
                       site == "WI" ~ 1,
                       site == "BA"  ~ 2)) %>% 
  ggplot(aes(x = jday, color = site, y = y)) +
  geom_point(aes(shape = site, fill = site),size = 4, alpha = 0.5) +
  scale_color_manual(values = colors[c(2,1,3,4)]) + 
  scale_fill_manual(values = colors[c(2,1,3,4)]) +
  ylab("") +
  geom_segment(aes(x = 96, xend = 189, y = 2, yend = 2), color = colors[1]) +
  geom_segment(aes(x = 87, xend = 198, y = 4, yend = 4), color = colors[2]) +
  geom_segment(aes(x = 73, xend = 208, y = 3, yend = 3), color = colors[4]) +
  geom_segment(aes(x = 55, xend = 154, y = 1, yend = 1), color = colors[3]) +
  geom_point(data = first_flower_site, aes(x = jday, y = c(4,2,1,3), fill = site),
             size = 4, shape = c(24,23,25,22), color = "black", stroke = 1) + 
  ylim(0.5,4.5) +
  annotate("text", x = 142.5, y = 3.4, label = "Cool seasonal (CH)", color = colors[4], size = 6) +
  annotate("text", x = 142.5, y = 4.4, label = "Cold aseasonal (SS)", color = colors[2], size = 6) +
  annotate("text", x = 140.5, y = 2.4, label = "Cool aseasonal (BA)", color = colors[1], size = 6) +
  annotate("text", x = 104.5, y = 1.4, label = "Hot seasonal (WI)", color = colors[3], size = 6) +
  theme_classic(base_size = 16) +
  scale_shape_manual(values = c(24,23,25,22)) +
  theme(legend.position = "none",
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank()) +
  xlab("Day of year") -> sampling_timeline 

png("figs/FigS1_SamplingFreq.png", height = 4, width = 6.6, units = "in", res = 300)
sampling_timeline
dev.off()
