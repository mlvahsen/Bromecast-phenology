# This file reads in and compiles the data needed for the flower timing analysis

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
## Create plant level dataset for flowering time ####
phen %>%
  filter(v %in% c("FG", "FP", "FB", "FX")) %>%
  group_by(plantID) %>%
  # Gets minimum day of flowering
  slice(which.min(jday)) -> phen_flower

# Update herbivory data at the level of the plant
phen %>% 
  filter(herbivory == "Y") %>% 
  pull(plantID) %>% unique() -> herbivory_plants 

phen_flower %>% 
  mutate(herbivory = ifelse(plantID %in% herbivory_plants, "Y", "N")) -> phen_flower

## Read in kinship data and subset flowering data to match ####
## Kinship stuff ####

# Read in data that matches kinship matrix position and genotype ID
kinshipIDs <- read_csv("~/Git/Bromecast/gardens/rawdata/cg_psuDTF.csv")

kinshipIDs %>% 
  # Two genotypes are currently absent from the kinship matrix
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
  mutate(genotype = as.factor(genotype))-> phen_flower_kin
# Right now this only drops 213 plants total

# And vice versa for kinship matrix
keeps <- which(rownames(kinship93BRTE) %in% unique(phen_flower_kin$genotype))

kinship93BRTE[keeps, keeps] -> kin

# Remove all other intermediate data sets
rm(list=setdiff(ls(), c("phen_flower_kin", "kin")))
