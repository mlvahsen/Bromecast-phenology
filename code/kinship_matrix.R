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

## Merge datasets together ####
phen <- rbind(phen_SS %>% dplyr::select(-tillers), phen_Boise)

## Calculate time to first flower ####

# Calculate when the earliest flowering day is for each plant. This procedure
# removes plants that do not flower of the course of the experiment.
phen %>%
  filter(v %in% c("FG", "FP", "FB")) %>%
  group_by(plantID) %>%
  # Gets minimum day of flowering
  slice(which.min(jday)) -> phen_flower

## Kinship stuff ####

# Read in data that matches kinship matrix position and genotype ID
kinshipIDs <- read_csv("~/Git/Bromecast/gardens/rawdata/cg_psuDTF.csv")

kinshipIDs %>% 
  filter(source != "Adler09" & source != "Shriver01") %>% 
  arrange(kinshipID)-> genotypes_93

# Read in kinship matrix
setwd("~/Git/Bromecast/gardens/rawdata/")
kinship93BRTE <- npyLoad("93BRTEcg.kinship.npy")
setwd("~/Git/Bromecast-phenology/")

# Put genotype numbers on rows and columns
colnames(kinship93BRTE) <- as.factor(genotypes_93$genotype)
rownames(kinship93BRTE) <- as.factor(genotypes_93$genotype)

# Subset phen data for only genotypes that we have kinship data for
phen_flower %>% 
  filter(genotype %in% genotypes_93$genotype) %>% 
  mutate(genotype = as.factor(genotype))-> phen_flower_sub

# And vice versa for kinship matrix
keeps <- which(rownames(kinship93BRTE) %in% unique(phen_flower_sub$genotype))

kinship93BRTE[keeps, keeps] -> kin

phen_flower %>% 
  filter(genotype %in% rownames(kin)) -> phen_flower_kin

## Fit linear model ####
start <- Sys.time()

gam_linear <- gam(jday ~ density*gravel*pc1 + density*gravel*pc2 +
                    pc1 + pc2 +
                    # Random intercept for block
                    s(block_unique, bs = 're') + 
                    # Random intercept for plot nested within block
                    s(plot_unique, bs = 're') +
                    s(genotype, bs = 're'),
                    #s(genotype, density, bs = 're') +
                    #s(genotype, site, bs = 're') +
                    #s(genotype, gravel, bs = "re"), method = "REML", 
                  data = phen_flowerSS)

end <- Sys.time()

phen_flower_kin %>% 
  filter(site == "SS") -> phen_flowerSS

brms_m1 <- brm(
  jday ~ 1 + density + gravel +
    pc1 + pc2 + (1 + density + gravel || gr(genotype, cov = Amat)) +
    (1 | block_unique) + (1 | plot_unique),
  data = phen_flowerBA,
  data2 = list(Amat = kin),
  family = gaussian(),
  chains = 2, cores = 1, iter = 1000
)

summary(brms_m1)

brms_m2 <- brm(
  jday ~ 1 + density*gravel +
    pc1 + pc2 + (1 + density + gravel || genotype) +
    (1 | block_unique) + (1 | plot_unique),
  data = phen_flower_kin,
  data2 = list(Amat = kin),
  family = gaussian(),
  chains = 1, cores = 1, iter = 1000
)

# Assess model fit
plot(brms_m2)
mcmc_plot(brms_m1, type = "acf")
summary(brms_m1)

# Calculate heritability
v_animal <- (VarCorr(brms_m1, summary = FALSE)$genotype$sd[,1])^2
v_r <- (VarCorr(brms_m1, summary = FALSE)$residual$sd)^2
h.bwt.1 <- as.mcmc(v_animal / (v_animal + v_r))
summary(h.bwt.1)
