# Calculate growing degree days (GDD) for each site

# Load libraries
library(tidyverse); library(prism); library(raster); library(geosphere); library(sf)

# Set graphics themes
theme_set(theme_bw(base_size = 16))

# Read GPS data for all genotypes
gps <- read_csv("~/Documents/Git/Bromecast Data/gardens/deriveddata/BioclimateOfOrigin_AllGenotypes.csv") %>% 
  dplyr::select(lon, lat, site_code) 

# Set download folder
prism_set_dl_dir("data/")

## Common garden site information ####

# Download mean temperature data for 2022
# get_prism_dailys(type = "tmean",
#                  minDate = "2021-10-01",
#                  maxDate = "2022-07-01",
#                  keepZip = F)

# Get tmax values for all 
to_slice <- prism_archive_subset("tmean", "daily", minDate="2021-10-01", maxDate = "2022-07-01")
stacked <- pd_stack(to_slice)
proj4string(stacked) <- CRS("+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")
df <- data.frame(rasterToPoints(stacked))
colnames(df) <- c("x", "y",1:274)

# Get gps of just cg sites
gps_sites <- gps %>% filter(site_code %in% c("SS", "BA", "WI", "CH"))

# Get closest prism point to each GPS point
store <- matrix(NA, nrow = nrow(gps_sites), ncol = 274)

for(i in 1:nrow(gps_sites)){
  out <- distm(gps_sites[i,c("lon", "lat")], df[,c("x", "y")], fun = distHaversine)
  store[i,] <- as.numeric(df[which.min(out),3:276])
}

# Match up prism data back to gps data frame
cbind(gps_sites, store) %>%
  gather(key = day, value = tmean, `1`:`274`) %>% 
  mutate(tmean = ifelse(tmean < 0, 0, tmean)) %>% 
  mutate(day = as.numeric(day)) -> site_temp

# Plot mean temp over the time period
# site_temp %>%
#   ggplot(aes(x = day, y = tmean, color = site_code)) +
#   geom_line()

# Create histograms of mean temp by site
# site_temp %>% 
#   mutate(tmean = case_when(tmean < 0 ~ 0,
#                            T ~ tmean)) %>%
#   ggplot(aes(x = tmean, fill = site_code)) + 
#   geom_density()+
#   facet_wrap(~site_code)

# Get average temperature values for each site
# site_temp %>% 
#   mutate(tmean = case_when(tmean < 0 ~ 0,
#                            T ~ tmean)) %>%  
#   group_by(site_code) %>% 
#   summarize(mean_temp = mean(tmean),
#             median_temp = median(tmean),
#             # Or calculate growing degree days (this creates the same
#             # relationships as means)
#             gdd = sum(tmean))

## Collection site information ####

# Download mean temperature data for 2022
# get_prism_normals(type="tmean",
#                   resolution = "4km",
#                   mon = c(1:6, 10:12),
#                   keepZip = FALSE)

`%notin%` <- Negate(`%in%`)

gps_collect <- gps %>% 
  filter(site_code %notin% c("SS", "BA", "CH", "WI"))

# Get tmax values for all 
to_slice <- prism_archive_subset("tmean", "monthly normals", resolution = "4km")
stacked <- pd_stack(to_slice)
proj4string(stacked) <- CRS("+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")
df <- data.frame(rasterToPoints(stacked))
colnames(df) <- c("x", "y",c(1:6, 10:12))

# Get closest prism point to each GPS point
store_collect <- matrix(NA, nrow = nrow(gps_collect), ncol = 9)

for(i in 1:nrow(gps_collect)){
  out <- distm(gps_collect[i,c("lon", "lat")], df[,c("x", "y")], fun = distHaversine)
  store_collect[i,] <- as.numeric(df[which.min(out),3:11])
}

# Match up prism data back to gps data frame
cbind(gps_collect, store_collect) %>% 
  gather(key = day, value = tmean, `1`:`9`) %>% 
  # Set any temps lower than 0 to be 0
  mutate(tmean = ifelse(tmean < 0, 0, tmean)) %>%
  # Drop observations in Canada
  filter(lat < 48) %>% 
  mutate(day = as.numeric(day)) -> collect_temp

# collect_temp %>% 
#   filter(lat < 48) %>% 
#   group_by(site_code) %>% 
#   summarize(mean_temp = mean(tmean),
#             median_temp = median(tmean),
#             # Or calculate growing degree days (this creates the same
#             # relationships as means)
#             gdd = sum(tmean)) %>% 
#   pull(mean_temp) %>% range()

# collect_temp %>% 
#   group_by(site_code) %>% 
#   summarize(mean_temp = mean(tmean),
#             median_temp = median(tmean),
#             # Or calculate growing degree days (this creates the same
#             # relationships as means)
#             gdd = sum(tmean)) %>% 
#   ggplot(aes(x = reorder(site_code, median_temp), y = median_temp)) +
#   geom_point()


