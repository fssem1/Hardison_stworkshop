library(tidyverse)
library(sf)
library(tinyVAST)
library(sdmTMB)
library(here)
library(DHARMa)
library(sp)
library(sdmTMBextra)
source(here::here("R/helpers.R"))
library(INLA)

# CRS
ak_crs <- "+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +datum=NAD83 +units=km +no_defs"

# when reading .kml files from google earth, there's
# a z attribute that should be dropped
pws_poly <- st_read(here::here("data/pws_poly.kml")) %>% 
  st_zm() %>% 
  st_transform(., ak_crs)

pws_poly_buff <- pws_poly %>% 
  # st_buffer(., dist = 10) %>%
  st_union() %>% 
  st_as_sf()

# Land polygon
ak_land <- rnaturalearthhires::states10 %>%
  filter(name == "Alaska") %>%
  st_union() %>% 
  st_transform(., ak_crs) %>% 
  st_intersection(., pws_poly) %>% 
  st_as_sf()


test <- 
  pws_poly_buff %>% 
  st_difference(., 
                ak_land) %>% 
  st_as_sf() %>% 
  mutate(test = "water")

## INLA mesh----
max.edge <- 20
offset <- 20
ak_land_sp <- test %>% as_Spatial()
custom_mesh <- inla.mesh.2d(boundary = ak_land_sp,
                            max.edge = c(max.edge, max.edge * 2),
                            offset = c(max.edge, offset))

custom_mesh$crs <- ak_crs
land.tri = inla.over_sp_mesh(ak_land_sp, y = custom_mesh, 
                             type = "centroid") 
num.tri = length(custom_mesh$graph$tv[, 1])
barrier.tri = setdiff(1:num.tri, land.tri)
poly.barrier = inla.barrier.polygon(custom_mesh, 
                                    barrier.triangles = barrier.tri)

plot(custom_mesh)
# generate observations in water

test_point <- test %>% 
  
  sf::st_sample(., 10000) %>%
  
  st_as_sf() %>% 
  
  sfc_as_cols(., names = c("longitude","latitude")) %>% 
  
  st_set_geometry(NULL)

# send to mesh

sdmTMB_mesh <- make_mesh(data = test_point,
                         
                         xy_cols = c("longitude","latitude"),
                         
                         cutoff = 0.5,
                         
                         range = 0.01)

# pass land polygon to add_barrier_mesh

sdmTMB_barrier_mesh <- 
  
  add_barrier_mesh(spde_obj = sdmTMB_mesh,
                   
                   barrier_sf = ak_land,
                   
                   plot = TRUE)

plot(sdmTMB_barrier_mesh)

# plot functions

mesh <- sdmTMB_barrier_mesh$mesh

tl <- length(mesh$graph$tv[, 1])

pos_tri <- matrix(0, tl, 2)

for (i in seq_len(tl)) {
  
  temp <- mesh$loc[mesh$graph$tv[i, ], ]
  
  pos_tri[i, ] <- colMeans(temp)[c(1, 2)]
  
}

mesh_sf <- as.data.frame(pos_tri)

mesh_sf$X <- mesh_sf$V1 * 1

mesh_sf$Y <- mesh_sf$V2 * 1

mesh_sf <- sf::st_as_sf(mesh_sf, crs = ak_crs, 
                        
                        coords = c("X", "Y"))

intersected <- sf::st_intersects(mesh_sf, ak_land)

water.triangles <- which(lengths(intersected) == 0)

land.triangles <- which(lengths(intersected) > 0)

water_verts <- 
  
  pos_tri[water.triangles, ] %>% 
  
  as.data.frame() %>% 
  
  st_as_sf(., coords = c("V1","V2"))

land_verts <- 
  
  pos_tri[land.triangles, ] %>% 
  
  as.data.frame() %>% 
  
  st_as_sf(., coords = c("V1","V2"))

plot(test)

plot(water_verts, add = TRUE)

plot(land_verts, add = TRUE, col = "red")


plot(poly.barrier, col = "#A0A2A1", add = T)
