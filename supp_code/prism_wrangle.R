# Calculate growing degree days (GDD) for each site

# Load libraries
library(tidyverse); library(prism); library(raster); library(geosphere); library(sf)

# Set graphics themes
theme_set(theme_bw(base_size = 16))

# Read GPS data for all genotypes
gps <- read_csv("~/Git/Bromecast/gardens/deriveddata/BioclimateOfOrigin_AllGenotypes.csv") %>% 
  dplyr::select(lon, lat, site_code) 

# Set download folder
prism_set_dl_dir("data/")

## Common garden site information ####

# Download mean temperature data for 2022
get_prism_dailys(type = "tmax",
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

# Plot mean temp over the time period
site_temp %>%
  ggplot(aes(x = day, y = tmean, color = site_code)) +
  geom_line() +
  scale_color_manual(values = colors[c(2,4,1,3)])

# Create histograms of mean temp by site
site_temp %>%
  mutate(tmean = case_when(tmean < 0 ~ 0,
                          T ~ tmean)) %>%
  ggplot(aes(x = tmean, fill = site_code, color = site_code)) +
  geom_density(alpha = 0, linewidth = 2) +
  scale_color_manual(values = colors[c(2,4,1,3)]) +
  scale_fill_manual(values = colors[c(2,4,1,3)])

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

write_csv(site_temp_summary, "~/Desktop/site_tmean.csv")


# Download mean temperature data for 2022
# get_prism_monthlys(type = "tmax",
#                  year = 2021,
#                  mon = 10:12,
#                  keepZip = F)
# 
# get_prism_monthlys(type = "tmax",
#                    year = 2022,
#                    mon = 1:6,
#                    keepZip = F)

# Get tmax values for all 
to_slice <- prism_archive_subset("tmean", "monthly")
stacked <- pd_stack(to_slice)
proj4string(stacked) <- CRS("+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")
df <- data.frame(rasterToPoints(stacked))
colnames(df) <- c("x", "y", "2021-10", "2021-11", "2021-12", "2022-01", "2022-02",
                  "2022-03", "2022-04", "2022-05", "2022-06")

# Get closest prism point to each GPS point
store <- matrix(NA, nrow = nrow(gps_sites), ncol = 9)

for(i in 1:nrow(gps_sites)){
  out <- distm(gps_sites[i,c("lon", "lat")], df[,c("x", "y")], fun = distHaversine)
  store[i,] <- as.numeric(df[which.min(out),3:11])
}

colnames(store) <- c("2021-10", "2021-11", "2021-12", "2022-01", "2022-02",
                  "2022-03", "2022-04", "2022-05", "2022-06")

# Match up prism data back to gps data frame
cbind(gps_sites, store) %>%
  gather(key = month, value = tmean, `2021-10`:`2022-06`)  -> site_temp

site_temp %>% 
  group_by(site_code) %>% 
  mutate(tmean = ifelse(tmean < 0, 0, tmean)) %>% 
  summarize(mean = mean(tmean))

library(ggforce)

site_temp %>% 
  ggplot(aes(x = month, y = tmean, color = site_code, group = site_code)) +
  ggforce::geom_link2(aes(group = site_code, color = after_stat(y > 0))) +
  geom_point(size = 3) +
  theme_classic(base_size = 16) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_color_manual(values = c(colors[2], colors[4], "gray87", colors[1], "black", colors[3])) +
  labs(x = "month in growing season", y = "mean monthly temperature (°C)") +
  geom_hline(aes(yintercept = 0), linetype = "dashed")



## Collection site information ####

# Download max temperature data for climate normals
get_prism_normals(type="tmax",
                  resolution = "4km",
                  mon = c(1:6),
                  keepZip = FALSE)

`%notin%` <- Negate(`%in%`)

gps_collect <- gps %>% 
  filter(site_code %notin% c("SS", "BA", "CH", "WI"))

# Get tmax values for all 
to_slice <- prism_archive_subset("tmean", "monthly normals", resolution = "4km")
stacked <- pd_stack(to_slice)
proj4string(stacked) <- CRS("+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")
df <- data.frame(rasterToPoints(stacked))
colnames(df) <- c("x", "y",c(1:6))

# Get closest prism point to each GPS point
store_collect <- matrix(NA, nrow = nrow(gps_collect), ncol = 6)

for(i in 1:nrow(gps_collect)){
  out <- distm(gps_collect[i,c("lon", "lat")], df[,c("x", "y")], fun = distHaversine)
  store_collect[i,] <- as.numeric(df[which.min(out),3:8])
}

# Match up prism data back to gps data frame
cbind(gps_collect, store_collect) %>% 
  gather(key = month, value = tmean, `1`:`6`) %>% 
  # Set any temps lower than 0 to be 0
  mutate(tmean = ifelse(tmean < 0, 0, tmean)) %>%
  # Drop observations in Canada
  filter(lat < 48) %>% 
  mutate(month = as.numeric(month)) -> collect_temp

collect_temp %>% 
  #filter(month < 5) %>% 
  group_by(lat, lon, site_code) %>% 
  summarize(tmax_mean = mean(tmean)) %>% 
  ggplot(aes(x = lon, y = lat, color = tmax_mean)) +
  geom_point(size = 3) +
  scale_color_distiller()

collect_temp %>% 
  group_by(lat, lon, site_code) %>% 
  summarize(tmean_mean = mean(tmean)) -> collect_temp_sum

write_csv(collect_temp_sum, "~/Desktop/genotype_tmean_norms.csv")

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


