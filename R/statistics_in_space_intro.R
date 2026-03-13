library(tidyverse)
library(sdmTMB)
library(sf)
library(patchwork)
library(gstat)
library(DHARMa)

# Load and process data------------------------------------------------------------

## Coordinate reference system for projecting CPUE data (Albers Equal Area)
## Ensure units are km
ak_crs <- "+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +datum=NAD83 +units=km +no_defs"

# bounding box for cropping Alaska land polygon
bb_poly <- st_as_sfc(st_bbox(c(
  xmin = -180,  xmax = -156,
  ymin =  50, ymax =  70), 
  crs = 4326))

## land polygon
ak_land <- rnaturalearthhires::states10 %>%
  filter(name == "Alaska") %>%
  st_union() %>% 
  st_intersection(bb_poly) %>% 
  st_transform(ak_crs)

## Helper functions
source(here::here("R/helpers.R"))

## load Red King Crab CPUE data
load(here::here("data/rkc_cpue.rdata"))

## process to sf and project
rkc_cpue_2010_sf <- rkc_cpue %>% filter(year == 2010,
                                     district == "BB") %>% 
  st_as_sf(., coords = c("longitude","latitude"), crs = 4326) %>% 
  st_transform(., ak_crs) 

## back to data frame
rkc_cpue_2010_df <- rkc_cpue_2010_sf %>% 
  sfc_as_cols(., names = c("longitude", "latitude")) %>% 
  st_set_geometry(NULL)

# map and CPUE obs------------------------------------------------------------------
p1.1 <- 
  ggplot(rkc_cpue_2010_sf) +
    geom_histogram(aes(cpue)) +
    scale_x_continuous(expand = c(0.01, 0.01)) +
    scale_y_continuous(expand = c(0.01, 0.01)) +
    theme_bw() +
    theme(axis.title = element_text(size = 14),
          title = element_text(size = 14)) +
    labs(title = "Bristol Bay red king crab CPUE histogram",
         x = "CPUE",
         y = "Count") 

p1.2 <- 
  ggplot() + 
    geom_sf(data = ak_land) +
    geom_sf(data = rkc_cpue_2010_sf) +
    coord_sf(ylim = c(550, 1100),
             xlim = c(-900, -175)) +
    theme_bw() +
    theme(axis.title = element_text(size = 14),
          title = element_text(size = 14)) +
    labs(title = "Bristol Bay survey stations",
         x = "Longitude",
         y = "Latitude") 

p1.2 + p1.1 + 
  plot_layout(nrow = 2,
              heights = c(1, 0.5))

# fit Tweedie intercept only GLM-----------------------------------------------------
mod1 <- sdmTMB(cpue ~ 1, 
               family = sdmTMB::tweedie(), 
               spatial = "off",
               data = rkc_cpue_2010_df)

summary(mod1)

# extract parameters + CIs
mod1_preds_df <- tidy(mod1) %>% 
  mutate(mu = exp(estimate),
         lwr = exp(conf.low),
         upr = exp(conf.high),
         Model = "GLM")

## map and CPUE obs with mean and CIs-------------------------------------------------
p1.3 <- 
  ggplot(rkc_cpue_2010_df) +
  geom_histogram(aes(cpue)) +
  geom_vline(data = mod1_preds_df,
             aes(xintercept = lwr),
             linetype = 2, color = "purple") +
  geom_vline(data = mod1_preds_df,
             aes(xintercept = upr),
             linetype = 2, color = "purple") +
  geom_vline(data = mod1_preds_df,
             aes(xintercept = mu), color = "purple") +
  scale_x_continuous(expand = c(0.01, 0.01)) +
  scale_y_continuous(expand = c(0.01, 0.01)) +
  theme_bw() +
  theme(axis.title = element_text(size = 14),
        title = element_text(size = 14)) +
  labs(title = "Bristol Bay red king crab CPUE histogram",
       x = "CPUE",
       y = "Count") 

p1.2 + p1.3 + 
  plot_layout(nrow = 2,
              heights = c(1, 0.5))

## evaluating residual independence---------------------------------------------------
resid_sim1 <- simulate(mod1, nsim = 500, type = "mle-mvn")
resid_dharm1 <- dharma_residuals(resid_sim1, mod1, return_DHARMa = TRUE)
DHARMa::testSpatialAutocorrelation(resid_dharm1, x = rkc_cpue_2010_df$longitude, y = rkc_cpue_2010_df$latitude)
rkc_cpue_2010_df$mod1_residuals <- resid_dharm1$scaledResiduals

## response residuals-----------------------------------------------------------------
rkc_cpue_2010_sf["response_residuals"] <- residuals(mod1, type = 'response')
resid_map_response <-
  ggplot() + 
  geom_sf(data = ak_land) +
  geom_sf(data = rkc_cpue_2010_sf, aes(color = response_residuals), size = 5) +
  scale_color_gradient2(low = scales::muted("red"),
                        high = scales::muted("blue"),
                        midpoint = 50) +
  coord_sf(ylim = c(550, 1100),
           xlim = c(-900, -175)) +
  theme_bw() +
  theme(axis.title = element_text(size = 14),
        title = element_text(size = 14),
        legend.position = "bottom") +
  labs(x = "Longitude",
       y = "Latitude",
       color = "Response\nresiduals") 
resid_map_response

## residual sample variogram----------------------------------------------------------
rkc_cpue_2010_spat <- rkc_cpue_2010_df 
sp::coordinates(rkc_cpue_2010_spat) <- c("longitude","latitude") # needs to be sp class for gstat
v1 <- gstat::variogram(mod1_residuals ~ 1, 
                       rkc_cpue_2010_spat)
sv_df <- 
  tibble(distance = v1$dist, 
       semivariance = v1$gamma,
       model = "GLM")

sv_plt <- 
  ggplot() + 
    geom_point(data = sv_df,
               aes(y = semivariance,
                   x = distance)) +
    geom_smooth(data = sv_df,
                aes(y = semivariance,
                    x = distance), se = F) +
    scale_x_continuous(expand = c(0.01, 0.01)) +
    scale_y_continuous(expand = c(0.01, 0.01)) +
    theme_bw() +
    theme(axis.title = element_text(size = 14),
          title = element_text(size = 14)) +
    labs(x = "Distance (km)",
         y = "Semivariance")

## map of scaled residuals----------------------------------------------------------
rkc_cpue_2010_sf["mod1_residuals"] <- resid_dharm1$scaledResiduals
resid_map <-
  ggplot() + 
  geom_sf(data = ak_land) +
  geom_sf(data = rkc_cpue_2010_sf, aes(color = mod1_residuals), size = 2) +
  scale_color_gradient2(low = scales::muted("red"),
                        high = scales::muted("blue"),
                        midpoint = 0.5) +
  coord_sf(ylim = c(550, 1100),
           xlim = c(-900, -175)) +
  theme_bw() +
  theme(axis.title = element_text(size = 14),
        title = element_text(size = 14),
        legend.position = "bottom") +
  labs(x = "Longitude",
       y = "Latitude",
       color = "Scaled\nresiduals") 


resid_map + sv_plt + 
  plot_layout(nrow = 1)

# fitting Tweedie GLMM with random intercept-----------------------------------------
## generate arbitrary grouping with kmeans
set.seed(100)
x <- rkc_cpue_2010_df %>% 
  dplyr::select(longitude, latitude) %>% 
  as.matrix()
cl <- kmeans(x, 15)
rkc_cpue_2010_df$group <- factor(cl$cluster)
rkc_cpue_2010_sf$group <- factor(cl$cluster)


## visualize clusters----------------------------------------------------------------
clust_map <- 
  ggplot() + 
  geom_sf(data = ak_land) +
  geom_sf(data = rkc_cpue_2010_sf, 
          aes(color= group), size = 4) +
  coord_sf(ylim = c(550, 1100),
           xlim = c(-900, -175)) +
  theme_bw() +
  theme(axis.title = element_text(size = 14),
        title = element_text(size = 14)) +
  labs(x = "Longitude",
       y = "Latitude") 
clust_map

## fit the GLMM----------------------------------------------------------------------
mod2 <- sdmTMB(cpue ~ 1 + (1|group), 
               family = sdmTMB::tweedie(), 
               spatial = "off",
               data = rkc_cpue_2010_df)

## DHARMa residuls------------------------------------------------------------------
resid_sim2 <- simulate(mod2, nsim = 500, type = "mle-mvn")
resid_dharm2 <- dharma_residuals(resid_sim2, mod2, return_DHARMa = TRUE)
DHARMa::testSpatialAutocorrelation(resid_dharm2, x = rkc_cpue_2010_df$longitude, y = rkc_cpue_2010_df$latitude)
rkc_cpue_2010_df$mod2_residuals <- resid_dharm2$scaledResiduals

## residual sample variogram--------------------------------------------------------
rkc_cpue_2010_spat <- rkc_cpue_2010_df
sp::coordinates(rkc_cpue_2010_spat) <- c("longitude","latitude")
v2 <- gstat::variogram(mod2_residuals ~ 1, 
                       rkc_cpue_2010_spat)
sv_df2 <- 
  tibble(distance = v2$dist, 
         semivariance = v2$gamma,
         model = "GLMM") %>% 
  bind_rows(., 
            sv_df)

sv_plt2 <- 
  ggplot() + 
  geom_point(data = sv_df2,
             aes(y = semivariance,
                 x = distance, color = model)) +
  geom_smooth(data = sv_df2,
              aes(y = semivariance,
                  x = distance, color = model), se = F) +
  scale_x_continuous(expand = c(0.01, 0.01)) +
  scale_y_continuous(expand = c(0.01, 0.01)) +
  theme_bw() +
  theme(axis.title = element_text(size = 14),
        title = element_text(size = 14),
        legend.position = "bottom") +
  labs(x = "Distance (km)",
       y = "Semivariance",
       color = "Model")

rkc_cpue_2010_sf["mod2_residuals"] <- resid_dharm2$scaledResiduals
resid_map2 <- 
  ggplot() + 
  geom_sf(data = ak_land) +
  geom_sf(data = rkc_cpue_2010_sf, aes(color = mod2_residuals), size = 2) +
  scale_color_gradient2(low = scales::muted("red"),
                        high = scales::muted("blue"),
                        midpoint = 0.5) +
  coord_sf(ylim = c(550, 1100),
           xlim = c(-900, -175)) +
  theme_bw() +
  theme(axis.title = element_text(size = 14),
        title = element_text(size = 14),
        legend.position = "bottom") +
  labs(x = "Longitude",
       y = "Latitude",
       color = "Scaled\nresiduals") 

resid_map2 + sv_plt2 + 
  plot_layout(nrow = 1)

## comparing estimates of RKC CPUE--------------------------------------------------
mod2_preds_df <- tidy(mod2) %>% 
  mutate(mu = exp(estimate),
         lwr = exp(conf.low),
         upr = exp(conf.high),
         Model = "GLMM") %>% 
  bind_rows(., mod1_preds_df)

## map and CPUE obs with mean and CIs-----------------------------------------------
glm_glmm_comp <-
  ggplot(rkc_cpue_2010_df) +
    geom_histogram(aes(cpue)) +
    geom_vline(data = mod2_preds_df,
               aes(xintercept = lwr,
                   color = Model),
               linetype = 2,linewidth = 1) +
    geom_vline(data = mod2_preds_df,
               aes(xintercept = upr,
                   color = Model),
               linetype = 2,linewidth = 1) +
    geom_vline(data = mod2_preds_df,
               aes(xintercept = mu,
                   color = Model),linewidth = 1) +
    scale_x_continuous(expand = c(0.01, 0.01)) +
    scale_y_continuous(expand = c(0.01, 0.01)) +
    theme_bw() +
    theme(axis.title = element_text(size = 14),
          title = element_text(size = 14)) +
    facet_wrap(~Model, nrow = 2)  +
    labs(x = "CPUE",
         y = "Count") 
glm_glmm_comp

# fitting Tweedie Spatial GLMM------------------------------------------------------
rkc_mesh <- make_mesh(data = rkc_cpue_2010_df, 
                      xy_cols = c("longitude", "latitude"), 
                      type = "cutoff", cutoff = 10)

mod3 <- sdmTMB(cpue ~ 1, 
               family = sdmTMB::tweedie(), 
               spatial = "on",
               mesh = rkc_mesh,
               data = rkc_cpue_2010_df)

## DHARMa residuals--------------------------
resid_sim3 <- simulate(mod3, nsim = 500, type = "mle-mvn")
resid_dharm3 <- dharma_residuals(resid_sim3, mod3, return_DHARMa = TRUE)
DHARMa::testSpatialAutocorrelation(resid_dharm3, x = rkc_cpue_2010_df$longitude, y = rkc_cpue_2010_df$latitude)
rkc_cpue_2010_df$mod3_residuals <- resid_dharm3$scaledResiduals

## residual sample variogram------------------------------------
rkc_cpue_2010_spat <- rkc_cpue_2010_df
sp::coordinates(rkc_cpue_2010_spat) <- c("longitude","latitude")
v3 <- gstat::variogram(mod3_residuals ~ 1, 
                       rkc_cpue_2010_spat)
sv_df3 <- 
  tibble(distance = v3$dist, 
         semivariance = v3$gamma,
         model = "Spatial GLMM") %>% 
  bind_rows(.,
            sv_df2)

rkc_cpue_2010_sf["mod3_residuals"] <- resid_dharm3$scaledResiduals
resid_map3 <- 
  ggplot() + 
  geom_sf(data = ak_land) +
  geom_sf(data = rkc_cpue_2010_sf, aes(color = mod3_residuals), size = 2) +
  scale_color_gradient2(low = scales::muted("red"),
                        high = scales::muted("blue"),
                        midpoint = 0.5) +
  coord_sf(ylim = c(550, 1100),
           xlim = c(-900, -175)) +
  theme_bw() +
  theme(axis.title = element_text(size = 14),
        title = element_text(size = 14),
        legend.position = "bottom") +
  labs(x = "Longitude",
       y = "Latitude",
       color = "Scaled\nresiduals") 

sv_plt3 <-
  ggplot() + 
  geom_point(data = sv_df3,
             aes(y = semivariance,
                 x = distance, color = model)) +
  geom_smooth(data = sv_df3,
              aes(y = semivariance,
                  x = distance, color = model), se = F) +
  scale_x_continuous(expand = c(0.01, 0.01)) +
  scale_y_continuous(expand = c(0.01, 0.01)) +
  theme_bw() +
  theme(axis.title = element_text(size = 14),
        title = element_text(size = 14),
        legend.position = "bottom") +
  labs(x = "Distance (km)",
       y = "Semivariance",
       color = "Model")
  
resid_map3 + sv_plt3 + 
  plot_layout(nrow = 1)

## Predict at unmeasured sampling locations-----------------------------------------
rkc_buff <- concaveman::concaveman(rkc_cpue_2010_sf, concavity = 3) %>% 
  st_buffer(dist = 10)

### gridded prediction region-------------------------------------------------------
rkc_gridded <- rkc_buff %>% 
  st_make_grid(., cellsize = 20) %>% 
  st_intersection(., rkc_buff) %>% 
  st_as_sf() %>% 
  mutate(id = 1:nrow(.))

### extract centroids and convert to df---------------------------------------------
rkc_pred_locs <- rkc_gridded %>% 
  st_centroid() %>% 
  st_as_sf() %>% 
  sfc_as_cols(., names = c("longitude", "latitude")) %>% 
  st_set_geometry(NULL)

### predict and map-------------------------------------------------------------------------
rkc_preds <- predict(mod3, newdata = rkc_pred_locs, type = "response")

### join predictions back with polygons
rkc_preds_sf <- rkc_gridded %>% 
  left_join(., rkc_preds)

pred_map <-
  ggplot() + 
    geom_sf(data = ak_land) +
    geom_sf(data = rkc_preds_sf, aes(fill = est, color = est)) +
    geom_sf(data = rkc_cpue_2010_sf %>% filter(cpue == 0), color = "red", alpha = 0.5) +
    geom_sf(data = rkc_cpue_2010_sf %>% filter(cpue > 0), aes(size = cpue)) +
    scale_color_gradient2(low = scales::muted("red"),
                          high = scales::muted("blue")) +
    scale_fill_gradient2(low = scales::muted("red"),
                          high = scales::muted("blue")) +
    coord_sf(ylim = c(550, 1100),
             xlim = c(-900, -175)) +
    guides(color = "none") +
    theme_bw() +
    theme(axis.title = element_text(size = 14),
          title = element_text(size = 14)) +
    labs(x = "Longitude",
         y = "Latitude",
         size = "Observed CPUE",
         fill = "Predicted CPUE") 

### predict on observation data
pred_on_obs <- predict(mod3, type = "response")
pred_comp_plt <- 
  ggplot() + 
    geom_point(data = pred_on_obs,
               aes(y = est, x = cpue)) +
    labs(y = "Predicted CPUE",
         x = "Observed CPUE") +
    theme_bw() + 
    theme(axis.title = element_text(size = 14),
          title = element_text(size = 14))

pred_map + pred_comp_plt +
  plot_layout(nrow = 1)
             
