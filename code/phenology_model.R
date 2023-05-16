# This code fits the phenology first to flower models to year 1 common garden
# data

## Load libraries ####
library(tidyverse); library(mgcv); library(gratia); library(geomtextpath);
library(here)
## Source code for genotype climate of origin ####
source(here("supp_code/climate_of_origin.R"))
## Read in all Boise data ####
# Read in derived phenology data
phen_Boise <- read_csv(url("https://raw.githubusercontent.com/pbadler/bromecast-data/main/gardens/deriveddata/Boise2022_growthphenology_with_harvest.csv?token=GHSAT0AAAAAACAWPIN25IM5XMROYGL6QXUKZDDX6XA"))
# Read in plant ID info
ids_Boise <- read_csv(url("https://raw.githubusercontent.com/pbadler/bromecast-data/main/gardens/deriveddata/Boise2022_plantID.csv?token=GHSAT0AAAAAACAWPIN2NZX6ATQMWZRHF7DWZDDX7KQ"))
# Read in garden treatment info
gardens <- read_csv(url("https://raw.githubusercontent.com/pbadler/bromecast-data/main/gardens/rawdata/garden_treatments.csv?token=GHSAT0AAAAAACAWPIN2IPEXXUZHZNZ5IQV6ZDDX75Q"))
# Read in flagging data
flags_Boise <- read_csv(url("https://raw.githubusercontent.com/pbadler/bromecast-data/main/gardens/deriveddata/Boise2022_flags.csv?token=GHSAT0AAAAAACAWPIN24VQIGPEOLS426FQOZDDYAOQ"))

# Merge together datasets
phen_id_Boise <- merge(phen_Boise, ids_Boise)

# Rename 'garden' column and remove cum_plot column
gardens %>% 
  mutate(site = garden) %>% 
  dplyr::select(-cum_plot, -garden) -> gardens_sub

# Merge together datasets
phen_id_garden_Boise <- merge(phen_id_Boise, gardens_sub)

# Merge together datasets
phen_Boise <- merge(phen_id_garden_Boise, flags_Boise)
phen_Boise <- merge(phen_Boise, genotype_PCclimate)

# Set appropriate factors for variables
phen_Boise %>% 
  mutate(block = as.factor(block),
         plot = as.factor(plot),
         growout = as.factor(growout),
         density = as.factor(density),
         gravel = as.factor(gravel),
         site = as.factor(site),
         genotype = as.factor(genotype)) %>% 
  mutate(plot_unique = as.factor(paste(site, block, plot, sep = "_")),
         block_unique = as.factor(paste(site, block, sep = "_")))-> phen_Boise

## Read in Sheep Station data ####

# Read in derived phenology data
phen_SS <- read_csv(url("https://raw.githubusercontent.com/pbadler/bromecast-data/main/gardens/deriveddata/SS2022_growthphenology_with_harvest.csv?token=GHSAT0AAAAAACAWPIN2VXFIJHDW2SRHMI5UZDDYA6Q"))
# Read in plant ID info
ids_SS <- read_csv(url("https://raw.githubusercontent.com/pbadler/bromecast-data/main/gardens/deriveddata/SS2022_plantID.csv?token=GHSAT0AAAAAACAWPIN25DABAIRUE4GB24NIZDDYBJA"))
# Read in flagging data
flags_SS <- read_csv(url("https://raw.githubusercontent.com/pbadler/bromecast-data/main/gardens/deriveddata/SS2022_flags.csv?token=GHSAT0AAAAAACAWPIN2JKAGEWJT6MUG7SJOZDDYCDA"))

# Merge together datasets
phen_id_SS <- merge(phen_SS, ids_SS)

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

## Bring together all data sets ####
# Remove tillers for SS
phen_SS %>% 
  dplyr::select(-tillers) -> phen_SS

phen <- rbind(phen_Boise, phen_SS)

# Remove other temporary objects
rm(list=setdiff(ls(), "phen"))

## Calculate time to first flower ####

# Calculate when the earliest flowering day is for each plant. This procedure
# removes plants that do not flower of the course of the experiment.
phen %>%
  filter(v %in% c("FG", "FP", "FB")) %>%
  group_by(plantID) %>%
  # Gets minimum day of flowering
  slice(which.min(jday)) -> phen_flower

# Get summary statistics of sampling effort
phen_flower %>%
  filter(site == "WI") %>% 
  group_by(jday, block) %>%
  summarize(n = n()) %>% 
  print(n = Inf)

phen_flower %>%
  filter(site == "SS") %>% 
  group_by(jday, block, plot) %>%
  summarize(n = n()) %>% 
  print(n = Inf)

# Reclassify sampling days to account for periods where not all blocks could be
# sampled on the same day
phen_flower %>% 
  mutate(jday = case_when(site == "WI" & jday == 131 ~ 130,
                          site == "WI" & jday == 153 ~ 151,
                          site == "SS" & jday %in% 124:125 ~ 123,
                          site == "SS" & jday %in% 139:140 ~ 138,
                          site == "SS" & jday %in% 152:154 ~ 151,
                          site == "SS" & jday == 166 ~ 165,
                          site == "SS" & jday == 187 ~ 186,
                          site == "SS" & jday == 189 ~ 188,
                          site == "SS" & jday %in% 193:195 ~ 192,
                          T ~ jday)) -> phen_flower

## Fit linear model ####
start <- Sys.time()

gam_linear <- gam(jday ~ density*gravel*pc1 + density*gravel*pc2 +
                    site*pc1 + site*pc2 +
                    # Random intercept for block
                    s(block_unique, bs = 're') + 
                    # Random intercept for plot nested within block
                    s(plot_unique, bs = 're') +
                    s(genotype, bs = 're') +
                    s(genotype, density, bs = 're') +
                    s(genotype, site, bs = 're') +
                    s(genotype, gravel, bs = "re"), method = "REML", 
                  data = phen_flower)

end <- Sys.time()

end - start
# Time difference of 6.9306 mins

# Look at model results
summary(gam_linear)

# Look at random effects
confint(gam_linear, parm = "s(block_unique)") %>% 
  ggplot(aes(x = reorder(block_unique, est), y = est)) +
  geom_point() +
  geom_segment(aes(y = lower, yend = upper, x = block_unique, xend = block_unique)) +
  xlab("genotype") +
  ylab("estimate") +
  geom_hline(aes(yintercept = 0), linetype = "dashed") +
  theme_bw(base_size = 14) +
  theme(axis.text.x = element_text(size = 8, angle = 45, hjust = 1))

# Show why (I think) block random intercept is getting estimated to be 0
phen_flower %>% 
  filter(site == "SS") %>% 
  group_by(block_unique) %>% 
  summarize(mean = mean(jday),
            sd = sd(jday),
            se = sd/sqrt(n())) -> summary_stats_blockWI

phen_flower %>% 
  filter(site == "SS") %>% 
  ggplot(aes(x = block_unique, y = jday)) +
  geom_count(aes(size = after_stat(prop), group = block_unique)) +
  scale_size_area(max_size = 10) +
  theme_bw(base_size = 14) +
  xlab("SS blocks") +
  geom_point(data = summary_stats_blockWI, aes(x = block_unique, y = mean),
             color = "blue", position = position_nudge(x = 0.2)) +
  geom_segment(data = summary_stats_blockWI, aes(x = block_unique,
                                                 xend = block_unique,
                                                 y = mean - 2*se, yend = mean + 2*se),
               position = position_nudge(x = 0.2), color = "blue") +
  labs(size = "proportion",
       y = "julian day")
  
## Fit ordinal model ####

# Create jday to ordered class data frame
tibble(jday = sort(unique(phen_flower$jday)),
       ord = 1:13) -> ords

# Merge with the rest of the data
merge(phen_flower, ords) -> phen_flower

gam_ord <- gam(ord ~ density*gravel*pc1 +
                density*gravel*pc2 + site*pc1 + site*pc2 +
                # Random intercept for block
                s(block_unique, bs = 're') + 
                # Random intercept for plot nested within block
                s(plot_unique, bs = 're') +
                s(genotype, bs = 're') +
                s(genotype, density, bs = 're') +
                s(genotype, site, bs = 're') +
                s(genotype, gravel, bs = "re"),
              family = ocat(R = 13),
              data = phen_flower)

variance_comp(gam_ord)

new_data <- expand.grid(density = unique(phen_flower$density),
                        plot_unique = unique(phen_flower$plot_unique),
                        genotype = unique(phen_flower$genotype),
                        gravel = unique(phen_flower$gravel)) %>% 
  separate(plot_unique, c("site", "block", "plot"), "_", remove = FALSE) %>% 
  mutate(block_unique = factor(paste(site, block, sep = "_")),
         site = as.factor(site),
         pc1 = mean(phen_flower$pc1),
         pc2 = mean(phen_flower$pc2)) %>% 
  dplyr::select(density, plot_unique, genotype, gravel, block_unique, site, pc1, pc2)

# Predict means across all possible combinations
predict(gam_ord, newdata = new_data, type = "response") -> m2_preds

# Add to exisiting data
cbind(new_data,
      p1 = m2_preds[,1],
      p2 = m2_preds[,2],
      p3 = m2_preds[,3],
      p4 = m2_preds[,4],
      p5 = m2_preds[,5],
      p6 = m2_preds[,6],
      p7 = m2_preds[,7],
      p8 = m2_preds[,8],
      p9 = m2_preds[,9],
      p10 = m2_preds[,10],
      p11 = m2_preds[,11],
      p12 = m2_preds[,12],
      p13 = m2_preds[,13]
      ) -> pred_data

as.data.frame(t(apply(m2_preds, 1, cumsum))) -> cum_preds
names(cum_preds) <- paste0("cp", 1:13, sep = "")

cbind(pred_data, cum_preds) -> pred_data

pred_data %>% 
  group_by(density, gravel, site) %>% 
  summarize(across(where(is.numeric), mean)) %>% 
  select(density, gravel, site, contains("cp")) %>% 
  gather(key = cat, value = cum_prob, cp1:cp13) %>% 
  mutate(jday = case_when(cat == "cp1" ~ 108,
                          cat == "cp2" ~ 123,
                          cat == "cp3" ~ 130,
                          cat == "cp4" ~ 138,
                          cat == "cp5" ~ 144,
                          cat == "cp6" ~ 147,
                          cat == "cp7" ~ 151,
                          cat == "cp8" ~ 154,
                          cat == "cp9" ~ 165,
                          cat == "cp10" ~ 182,
                          cat == "cp11" ~ 186,
                          cat == "cp12" ~ 188,
                          cat == "cp13" ~ 192),
         treatment = paste(density, gravel, sep = "_"),
         density = ifelse(density == "hi", "high", "low")) %>% 
  ggplot(aes(x = jday, y = cum_prob, group = treatment, color = gravel, shape = density, linetype = density)) +
  geom_line() +
  geom_point(size = 3) +
  ylab("cumulative probability of flowering by date") +
  theme_classic(base_size = 14) +
  xlab("julian day") +
  scale_color_manual(values = c("black", "gray70")) +
  scale_x_continuous(breaks = c(108,123,130,138,144,147,151,154,165,182,186,188,192)) +
  facet_wrap(~site) +
  theme(axis.text.x = element_text(size = 8, angle = 45, hjust = 1))

phen_flower %>% 
  filter(site == "WI") %>% 
  group_by(plot_unique, gravel, density) %>% 
  summarize(n = n())



