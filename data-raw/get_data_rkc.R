# packages
library(tidyverse)
library(crabpack) #devtools::install_github("AFSC-Shellfish-Assessment-Program/crabpack")

# download and process NMFS EBS BTS data and calculate CPUE using crabpack
species <- "RKC"
region1 <- "EBS"
years <- 2005:2024
channel <- "API"

# queries----
# query works, but breaks calc_cpue if region > 1
specimen_data_ebs_rkc <- crabpack::get_specimen_data(species = species,
                                                 region = region1, # Eastern Bering Sea and Northern Bering Sea
                                                 years = years,
                                                 channel = channel)

# get haul data for (gear) temperature and depth
haul_df_rkc <- specimen_data_ebs_rkc$haul %>% 
  rename_all(str_to_lower) %>% 
  dplyr::rename(latitude = mid_latitude,
                longitude = mid_longitude) %>%
  dplyr::select(bottom_depth,
                temp = gear_temperature,
                longitude,
                latitude,
                year,
                region)

# calculate CPUE
ebs_cpue_rkc <- 
  calc_cpue(crab_data = specimen_data_ebs_rkc,
            species = species,
            region = region1,
            year = years,
            crab_category = "legal_male")

# combine CPUE
rkc_cpue <- ebs_cpue_rkc %>% 
  rename_all(str_to_lower) %>% 
  dplyr::select(year, region, station_id, latitude, longitude, category, district, cpue) %>% 
  mutate(cpue = cpue/3.429) %>% # convert from catch (numbers) per sq. nm to catch per sq. km
  left_join(.,
            haul_df_rkc)

save(rkc_cpue,
     file = here::here("data/rkc_cpue.rdata"))
