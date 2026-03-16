library(tidyverse)
library(sf)
library(sdmTMB)
library(here)
library(INLA) #
library(gstat) #
library(mgcv) 
library(sp)
library(DHARMa) #
library(sdmTMBextra) #remotes::install_github("pbs-assess/sdmTMBextra", dependencies = TRUE)

# load processed data from R/data_processing_and_visualization
load(here::here("data/processed_hair_crab_data.rdata"))

# load hair crab data (includes design-based index)
load(here::here("data/hair_cpue.rdata"))

# source helper functions
source(here::here("R/helpers.R"))

# Convert sf to data.frame for use with sdmTMB
hair_cpue <- 
  hair_cpue_sf_filt %>% 
  sfc_as_cols(., names = c("longitude", "latitude")) %>% # turn points into lat/lon columns
  st_set_geometry(NULL) %>% #convert to data.frame
  filter(year >= 2014) %>% 
  mutate(district = factor(district))


# Approaches within sdmTMB vignette------
# Extremely simple cutoff:
sdmTMB_mesh1 <- make_mesh(hair_cpue, c("longitude", "latitude"), 
                          cutoff = 5, type = "cutoff")
plot(sdmTMB_mesh1)

# Using a k-means algorithm to assign vertices:
sdmTMB_mesh2 <- make_mesh(hair_cpue, c("longitude", "latitude"), 
                          n_knots = 150, type = "kmeans")
plot(sdmTMB_mesh2)


# CRS
ak_crs <- "+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +datum=NAD83 +units=km +no_defs"

# add in 2020 data
m1.2 <- sdmTMB::sdmTMB(
    formula = cpue ~ 1,
    extra_time = c(2020),
    spatial = "on",
    spatiotemporal = "ar1",
    time = "year",
    anisotropy = TRUE,
    mesh = sdmTMB_mesh2,  # simple mesh for speed
    family = tweedie(link = "log"),
    data = hair_cpue)

m1.2

# not including 2020
m1 <- sdmTMB::sdmTMB(
  formula = cpue ~ 0 + as.factor(year),
  spatial = "on",
  spatiotemporal = "iid",
  time = "year",
  mesh = sdmTMB_mesh2, # simple mesh for speed
  family = tweedie(link = "log"),
  data = hair_cpue)
save(m1, file = here::here("data/m1.rdata"))

m1


pred_out_ar <- predict(m1.2.ar, newdata = pred_df, type = "link",
                                return_tmb_object = TRUE) 
pred_out_iid <- predict(m1.2, newdata = pred_df %>% filter(year!=2020), type = "link",
                                return_tmb_object = TRUE) 