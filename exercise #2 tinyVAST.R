library(tidyverse)
library(sf)
library(tinyVAST)
library(sdmTMB)
library(here)
library(DHARMa)

# load processed data from R/data_processing_and_visualization
load(here::here("data/processed_hair_crab_data.rdata"))

# load hair crab data
load(here::here("data/hair_cpue.rdata"))

# source helper functions
source(here::here("R/helpers.R"))

# Convert sf to data.frame for use with sdmTMB
hair_cpue <- 
  hair_cpue_sf_filt %>% 
  sfc_as_cols(., names = c("longitude", "latitude")) %>% # turn points into lat/lon columns
  st_set_geometry(NULL) %>% #convert to data.frame
  filter(year >= 2014) %>% 
  mutate(district = ifelse(district == "N. EBS/NBS", "NEBS_NBS", district),
         var = "hair_crab",
  ) %>% 
  as.data.frame()

# CRS
ak_crs <- "+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +datum=NAD83 +units=km +no_defs"

# load fitted sdmTMB model
load(here::here("data/m1_2.rdata"))

# We will use the same mesh for tinyVAST
sdmTMB_mesh2 <- make_mesh(hair_cpue, c("longitude", "latitude"), 
                          n_knots = 150, type = "kmeans")
tinyVAST_mesh <- sdmTMB_mesh2$mesh

space_term <- "
  hair_crab <-> hair_crab, sd_space
"

spacetime_term_ar1 <- "
  hair_crab -> hair_crab, 1, rho 
  hair_crab <-> hair_crab, 0, sd_spacetime
"

# AR1------
m1_tv_ar1 <- tinyVAST(
  space_term = space_term,
  spacetime_term = spacetime_term_ar1,
  data = hair_cpue,
  formula = cpue ~ 1,
  spatial_domain = tinyVAST_mesh,
  family = tweedie(),
  time_column = "year",
  space_columns = c("longitude","latitude"),
  variable_column = "var",
  control = tinyVASTcontrol(
    use_anisotropy = TRUE
  )
)
m1_tv_ar1

# prediction grid
# polygon for the EBS
cellsize <- 20 #km x 20 km

# polygon for the EBS
ebs_domain <- 
  st_read(here::here("data/SAP_layers.gdb"),
          layer = "EBS_grid",
          quiet = TRUE) %>% 
  st_transform(., ak_crs) %>% 
  st_union()

# re-grid it to our spatial resolution
ebs_grid <- ebs_domain %>% 
  st_make_grid(., cellsize  = cellsize) %>% 
  st_intersection(., hair_loc_buffers) %>%
  st_as_sf() %>% 
  mutate(id = 1:nrow(.),
         area = as.numeric(st_area(.)))

# area vector to pass to integrate_output()
area_vec <- ebs_grid$area

# expand over time
pred_grid <- ebs_grid %>% 
  tidyr::expand_grid(.,
                     year = 2014:2030) %>% 
  st_as_sf() %>% 
  filter(year != 2020) # COVID year - no survey so should not be in final data

# final outputs
year_seq <- 2025:2030
pred_df <- pred_grid %>% 
  st_centroid() %>% 
  sfc_as_cols(., names = c("longitude","latitude")) %>% 
  st_set_geometry(NULL) %>% 
  as.data.frame()%>% 
  filter(year %in% year_seq)

pred_df$projection <-
tinyVAST::project(object = m1_tv_ar1,
                  extra_times = year_seq,
                  newdata = pred_df)



