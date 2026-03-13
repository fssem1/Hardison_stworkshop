# packages
library(tidyverse)
library(crabpack) #devtools::install_github("AFSC-Shellfish-Assessment-Program/crabpack")

# download and process NMFS EBS BTS data and calculate CPUE using crabpack
species <- "HAIR"
region1 <- "EBS"
region2 <- "NBS"
years <- 2005:2024
channel <- "API"

# queries----
# query works, but breaks calc_cpue if region > 1
specimen_data_ebs <- crabpack::get_specimen_data(species = species,
                                             region = region1, # Eastern Bering Sea and Northern Bering Sea
                                             years = years,
                                             channel = channel)

# calculate design-based index for legal male hair crab
hc_design_index_ebs <- crabpack::calc_bioabund(crab_data = specimen_data_ebs,
                                           species = species,
                                           region = region1,
                                           spatial_level = "region",
                                           crab_category = "legal_male")

specimen_data_nbs <- crabpack::get_specimen_data(species = species,
                                                 region = region2, # Eastern Bering Sea and Northern Bering Sea
                                                 years = years,
                                                 channel = channel)

# get haul data for (gear) temperature and depth
haul_ebs <- specimen_data_ebs$haul %>% 
  rename_all(str_to_lower) %>% 
  dplyr::rename(latitude = mid_latitude,
                longitude = mid_longitude) %>%
  dplyr::select(bottom_depth,
                temp = gear_temperature,
                longitude,
                latitude,
                year,
                region)


haul_nbs <- specimen_data_nbs$haul %>% 
  rename_all(str_to_lower) %>% 
  dplyr::rename(latitude = mid_latitude,
                longitude = mid_longitude) %>%
  dplyr::select(bottom_depth,
                temp = gear_temperature,
                longitude,
                latitude,
                year,
                region)

# combine
haul_df <- bind_rows(haul_ebs,
                     haul_nbs)

# calculate CPUE
ebs_cpue <- 
  calc_cpue(crab_data = specimen_data_ebs,
           species = species,
           region = region1,
           year = years,
           crab_category = "legal_male")

nbs_cpue <- 
  calc_cpue(crab_data = specimen_data_nbs,
            species = species,
            region = region2,
            year = years,
            crab_category = "legal_male")

# combine CPUE
hair_cpue <- bind_rows(ebs_cpue,
                       nbs_cpue) %>% 
  rename_all(str_to_lower) %>% 
  dplyr::select(year, region, station_id, latitude, longitude, category, district, cpue) %>% 
  mutate(cpue = cpue/3.429) %>% # convert from catch (numbers) per sq. nm to catch per sq. km
  left_join(.,
            haul_df)

# calculate CPUE
ebs_cpue_female <- 
  calc_cpue(crab_data = specimen_data_ebs,
            species = species,
            region = region1,
            year = years,
            crab_category = "mature_female")

nbs_cpue_female <- 
  calc_cpue(crab_data = specimen_data_nbs,
            species = species,
            region = region2,
            year = years,
            crab_category = "mature_female")

# combine CPUE
hair_cpue_female <- bind_rows(ebs_cpue_female,
                       nbs_cpue_female) %>% 
  rename_all(str_to_lower) %>% 
  dplyr::select(year, region, station_id, latitude, longitude, category, district, cpue) %>% 
  mutate(cpue = cpue/3.429) %>% # convert from catch (numbers) per sq. nm to catch per sq. km
  left_join(.,
            haul_df)
  
save(hair_cpue,
     hair_cpue_female,
     hc_design_index_ebs,
     file = here::here("data/hair_cpue.rdata"))

# Notes from EBS report (Zacher et al 2024)

# In this report, legal male hair crab (Erimacrus isenbeckii) are defined as > 3.25 inches CW (≥ 83 mm CL), which was specified in the previous 
# Pribilof District fishery; the female hair crab biomass estimate is presented for all sizes and maturity states combined. Hair crab were 
# caught at 46 of the 349 stations throughout all districts combined on the survey (Fig. 105). The 2024 biomass estimate of legal males was 
# 923 ± 358 t (1.7 ± 0.7 million crab) and 579 ± 250 t (1.8 ± 0.9 million crab) for sub-legal males (Tables 36 and 37; Fig. 101). Male hair 
# crab primarily occurred along the 50 m isobath and into Bristol Bay (Figs. 102, 103, and 105). The female hair crab biomass estimate was 
# 454 ± 274 t (1.3 ± 0.7 million crab; Tables 36 and 37; Fig. 101). Females were primarily distributed near the 50 m isobath, in Bristol Bay,
# and around Saint Paul Island (Figs. 104 and 105).
