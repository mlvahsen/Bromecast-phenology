# Calculate growing degree days (GDD) for each site

# Load libraries
library(tidyverse); library(prism); library(raster); library(geosphere)

# Read GPS data for all genotypes
gps <- read_csv("~/Git/Bromecast/gardens/deriveddata/BioclimateOfOrigin_AllGenotypes.csv") %>% 
  dplyr::select(lon, lat, site_code) %>% 
  filter(site_code %in% c("SS", "BA", "WI", "CH"))

# Set download folder
prism_set_dl_dir("data/")

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
colnames(df) <- c("x", "y",1:244)

# Get closest prism point to each GPS point
store <- matrix(NA, nrow = nrow(gps), ncol = 244)

for(i in 1:nrow(gps)){
  out <- distm(gps[i,c("lon", "lat")], df[,c("x", "y")], fun = distHaversine)
  store[i,] <- as.numeric(df[which.min(out),3:246])
}

