# Get PRISM data for each site (across the growing season) and genotype (across 30
# year climate norm)

# Load libraries
library(tidyverse); library(prism); library(raster); library(geosphere); library(sf)

# Set graphics themes
theme_set(theme_bw(base_size = 16))

# Read GPS data for all genotypes from main Bromecast repository
gps <- read_csv("https://raw.githubusercontent.com/pbadler/bromecast-data/main/gardens/deriveddata/BioclimateOfOrigin_AllGenotypes.csv") %>% 
  dplyr::select(lon, lat, site_code) 

# Set download folder for PRISM data
prism_set_dl_dir("data/")

## Common garden site information ####

# Download daily mean temperature data for 2022
get_prism_dailys(type = "tmean",
                 minDate = "2021-10-01",
                 maxDate = "2022-07-01",
                 keepZip = F)

# Get tmean values for all 
to_slice <- prism_archive_subset("tmean", "daily", minDate="2021-10-01", maxDate = "2022-06-30")
stacked <- pd_stack(to_slice)
proj4string(stacked) <- CRS("+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")
df <- data.frame(rasterToPoints(stacked))
colnames(df) <- c("x", "y",1:273)

# Get gps of just cg sites
gps_sites <- gps %>% filter(site_code %in% c("SS", "BA", "WI", "CH"))

# Get closest prism point to each GPS point
store <- matrix(NA, nrow = nrow(gps_sites), ncol = 273)

for(i in 1:nrow(gps_sites)){
  out <- distm(gps_sites[i,c("lon", "lat")], df[,c("x", "y")], fun = distHaversine)
  store[i,] <- as.numeric(df[which.min(out),3:275])
}

# Match up prism data back to gps data frame
cbind(gps_sites, store) %>%
  gather(key = day, value = tmean, `1`:`273`) %>% 
  mutate(tmean = ifelse(tmean < 0, 0, tmean)) %>% 
  mutate(day = as.numeric(day)) -> site_temp

# Set site colors
colors <- c("#88CCEE", "#AA4499", "#DDCC77", "#44AA99")

# Get average temperature values for each site
site_temp %>%
  mutate(tmean = case_when(tmean < 0 ~ 0,
                           T ~ tmean)) %>%
  group_by(site_code) %>%
  summarize(mean_temp = mean(tmean),
            median_temp = median(tmean),
            # Or calculate growing degree days (this creates the same
            # relationships as means)
            gdd = sum(tmean)) -> site_temp_summary

write_csv(site_temp_summary, "supp_data/site_tmean.csv")

## Collection site information - tmax ####

# Download max temperature data for climate normals
get_prism_normals(type="tmax",
                  resolution = "4km",
                  mon = 1:4,
                  keepZip = FALSE)

# Create %notin% operator
`%notin%` <- Negate(`%in%`)

# Get just gps collection sites
gps_collect <- gps %>% 
  filter(site_code %notin% c("SS", "BA", "CH", "WI"))

# Get tmax values for all 
to_slice <- prism_archive_subset("tmax", "monthly normals", resolution = "4km")
stacked <- pd_stack(to_slice)
proj4string(stacked) <- CRS("+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")
df <- data.frame(rasterToPoints(stacked))
colnames(df) <- c("x", "y", 1:4)

# Get closest prism point to each GPS point
store_collect <- matrix(NA, nrow = nrow(gps_collect), ncol = 4)

for(i in 1:nrow(gps_collect)){
  out <- distm(gps_collect[i,c("lon", "lat")], df[,c("x", "y")], fun = distHaversine)
  store_collect[i,] <- as.numeric(df[which.min(out),3:6])
}

# Match up prism data back to gps data frame
cbind(gps_collect, store_collect) %>% 
  gather(key = month, value = tmax, `1`:`4`) %>% 
  # Set any temps lower than 0 to be 0
  mutate(tmax = ifelse(tmax < 0, 0, tmax)) %>%
  # Drop observations in Canada because PRISM data doesn't work
  filter(lat < 48) %>% 
  mutate(month = as.numeric(month)) -> collect_temp

collect_temp %>% 
  group_by(lat, lon, site_code) %>% 
  summarize(tmax_mean = mean(tmax)) -> collect_temp_sum

# Save PRISM data for genotypes
write_csv(collect_temp_sum, "supp_data/genotype_tmax_norms.csv")

## Collection site information - tmean ####

# Download max temperature data for climate normals
get_prism_normals(type="tmean",
                  resolution = "4km",
                  mon = c(1:6,10:12),
                  keepZip = FALSE)

# Get tmean values for all 
to_slice <- prism_archive_subset("tmean", "monthly normals", resolution = "4km")
stacked <- pd_stack(to_slice)
proj4string(stacked) <- CRS("+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")
df <- data.frame(rasterToPoints(stacked))
colnames(df) <- c("x", "y", c(1:6,10:12))

# Get closest prism point to each GPS point
store_collect <- matrix(NA, nrow = nrow(gps_collect), ncol = 9)

for(i in 1:nrow(gps_collect)){
  out <- distm(gps_collect[i,c("lon", "lat")], df[,c("x", "y")], fun = distHaversine)
  store_collect[i,] <- as.numeric(df[which.min(out),3:11])
}

# Match up prism data back to gps data frame
cbind(gps_collect, store_collect) %>% 
  gather(key = month, value = tmean, `1`:`9`) %>% 
  # Set any temps lower than 0 to be 0
  mutate(tmean = ifelse(tmean < 0, 0, tmean)) %>%
  # Drop observations in Canada because PRISM data doesn't work
  filter(lat < 48) %>% 
  mutate(month = as.numeric(month)) -> collect_temp

collect_temp %>% 
  group_by(lat, lon, site_code) %>% 
  summarize(tmean_mean = mean(tmean)) -> collect_temp_sum

# Save PRISM data for genotypes
write_csv(collect_temp_sum, "supp_data/genotype_tmean_norms.csv")
