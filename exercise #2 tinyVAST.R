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

start_year <- 2010

end_year <- 2023

# combine male and female data and convert to sf

hair_cpue_sf_mf <- 
  
  hair_cpue_sf_all %>% 
  
  bind_rows(., hair_cpue_female %>% 
              
              st_as_sf(., coords = c("longitude", "latitude"), 
                       
                       crs = 4326) %>% 
              
              st_transform(ak_crs)) %>% 
  
  dplyr::rename(var = category) %>% 
  
  filter(between(year, start_year, end_year))

# create buffer around positive catches to drop false zeros before fitting model

hair_loc_buffers_mf <- hair_cpue_sf_mf %>% # across all trawls
  
  filter(cpue > 0) %>% 
  
  concaveman::concaveman(., concavity = 3) %>% #Draw concave hull around positive catches
  
  st_buffer(., dist = 45) %>%  # extend buffer
  
  st_difference(ak_land) #remove sections of polygon intersecting land

# intersect CPUE observations with buffer object

hair_cpue_mf <- hair_cpue_sf_mf %>%
  
  mutate(district = ifelse(district %in% c("ALL", "NORTH"), 
                           
                           "N. EBS/NBS",
                           
                           district)) %>% 
  
  st_intersection(., hair_loc_buffers_mf) %>% 
  
  sfc_as_cols(., names = c("longitude", "latitude")) %>% 
  
  st_set_geometry(NULL) %>% 
  
  as.data.frame()

# define new mesh

tinyVAST_mesh2 <- sdmTMB::make_mesh(hair_cpue_mf,
                                    
                                    xy_cols = c("longitude","latitude"),
                                    
                                    type = "kmeans",
                                    
                                    n_knots = 150)$mesh

spacetime_term_vast <- "

  # autoregressive effects

  legal_male -> legal_male, 1, b11

  mature_female -> mature_female, 1, b22

  # cross lagged effects

  mature_female -> legal_male, 1, b21

  legal_male -> mature_female, 1, b12

  # spatiotemporal variances (seperate)

  legal_male <-> legal_male, 0, sd_spacetime1

  mature_female <-> mature_female, 0, sd_spacetime2

"

# fit model

m2_tv <- tinyVAST(
  
  spacetime_term = spacetime_term_vast,
  
  data = hair_cpue_mf,
  
  time_column = "year",
  
  space_columns = c("longitude","latitude"),
  
  variable_column = "var",
  
  formula = cpue ~ 0 + var,
  
  spatial_domain = tinyVAST_mesh2,
  
  family = tweedie() )


spacetime_term_vast2 <- "

  # autoregressive effects

  legal_male -> legal_male, 1, b11

  # spatiotemporal variances (seperate)

  legal_male <-> legal_male, 0, sd_spacetime1

"

hair_cpue_male <- hair_cpue_mf %>% filter(var == "legal_male")

m3_tv <- tinyVAST(
  
  spacetime_term = spacetime_term_vast2,
  
  data = hair_cpue_male,
  
  time_column = "year",
  
  space_columns = c("longitude","latitude"),
  
  variable_column = "var",
  
  formula = cpue ~ 1,
  
  spatial_domain = tinyVAST_mesh2,
  
  family = tweedie() ,
  
  control = tinyVASTcontrol(trace = getOption("tinyVAST.trace", 1),
                            
                            silent = getOption("tinyVAST.silent", F)))

# get M/F CPUE data for 2024

hair_cpue_mf_2024 <- 
  
  hair_cpue_sf_all %>% 
  
  bind_rows(., hair_cpue_female %>% 
              
              st_as_sf(., coords = c("longitude", 
                                     
                                     "latitude"), crs = 4326) %>% 
              
              st_transform(ak_crs)) %>% 
  
  dplyr::rename(var = category) %>% 
  
  filter(year == 2024) %>% 
  
  sfc_as_cols(., names = c("longitude","latitude")) %>% 
  
  st_set_geometry(NULL) %>% 
  
  as.data.frame()

# get male only CPUE data for 2024

hair_cpue_m_2024 <- hair_cpue_mf_2024 %>% 
  
  filter(var == "legal_male")

extra_time <- 2024

hair_cpue_mf_2024$projection <- project(m2_tv, 
                                        
                                        newdata = hair_cpue_mf_2024,
                                        
                                        extra_times = extra_time)

hair_cpue_m_2024$projection <- project(m3_tv, 
                                       
                                       newdata = hair_cpue_m_2024,
                                       
                                       extra_times = extra_time)
