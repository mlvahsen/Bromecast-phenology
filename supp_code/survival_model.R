# Survival models

# This code fits th


## Calculate time to first flower ####

# Figure out survival model
phen %>%
  group_by(plantID) %>%
  slice(which.max(jday)) %>%
  mutate(survived = ifelse(live == "Y", 1, 0)) %>%
  dplyr::select(plantID, site, block, block_unique, plot, plot_unique, frostheave_date, density, gravel, genotype, survived, pc1, pc2) %>%
  distinct() -> phen_survival

mod_survival <- glmer(survived ~ density*gravel*pc1 + density*gravel*pc2 +
                        site*pc1 + site*pc2 + (1|block_unique) + (1|plot_unique)+(1+density+gravel+site|genotype),
                      family = "binomial", data = phen_survival)

# Figure out flower binary model
phen_survival %>%
  filter(survived == 1) %>%
  pull(plantID) -> survived_plants

phen %>%
  filter(v %in% c("FG", "FB", "FP", "FX") & plantID %in% survived_plants) %>%
  group_by(plantID) %>%
  slice(which.min(jday)) %>%
  ungroup() %>%
  mutate(flowered = 1) -> flowered_plants

# Subset to get plants that did not flower
`%notin%` <- Negate(`%in%`) # Not in operator

phen %>%
  filter(plantID %in% survived_plants) %>%
  filter(plantID %notin% flowered_plants$plantID) %>%
  group_by(plantID) %>%
  slice(which.min(jday)) %>%
  ungroup() %>%
  mutate(flowered = 0) -> unflowered_plants

# Join together flowered and unflowered data sets
phen_flowering <- rbind(flowered_plants, unflowered_plants)

mod_flowered <- glmer(flowered ~ density*gravel*pc1 + density*gravel*pc2 +
                        site*pc1 + site*pc2 + (1|block_unique) + (1|plot_unique)+(1+density+gravel+site|genotype),
                      family = "binomial", data = phen_flowering)