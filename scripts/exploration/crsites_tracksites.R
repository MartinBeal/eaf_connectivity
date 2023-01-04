### Color ringing sites - derive from data ---------------------------------

pacman::p_unload(plyr)
pacman::p_load(dplyr, magrittr, stringr, sf, ggplot2)

crdat <- readRDS("data/analysis/ringing/cr_merge.rds") # from cr-schemes

# vs <- filter(crdat, scheme == "germ_VSalewski")

nas <- filter(crdat, is.na(obstype))
table(nas$scheme)

firsts <- crdat %>%
  arrange(bird_id, date) %>% group_by(bird_id) %>%
  slice(1) %>% 
  filter(obstype %in% "N" | is.na(obstype))

## 
firsts %>%
  sf::st_as_sf(coords = c("longitude", "latitude"),
               crs = 4326, agr = "constant") %>%
  mapview::mapview()


### Tracking deployment sites ---------------------------------

pacman::p_unload(plyr)
pacman::p_load(dplyr, magrittr, stringr, sf, ggplot2)

## tracking data
TD <- readRDS("data/analysis/tracking/PTT_GPS_mconn_12h.rds")

firsts <- TD %>%
  arrange(id, timestamp) %>% 
  group_by(id) %>%
  filter(argos_lc %in% c("G", "1", "2", "3", "g")) %>% 
  slice(1)


## 
firsts %>%
  sf::st_as_sf(coords = c("longitude", "latitude"),
               crs = 4326, agr = "constant") %>%
  mapview::mapview()

##
firsts_bin <- firsts %>% 
  mutate(
    latitude  = round(latitude, 1),
    longitude = round(longitude, 1)
  ) %>% group_by(latitude, longitude) %>% summarise(
    n_id = n_distinct(id)
  )


## 
firsts_bin %>%
  sf::st_as_sf(coords = c("longitude", "latitude"),
               crs = 4326, agr = "constant") %>%
  mapview::mapview()


## Static map: project for prettier map --------------------------------------
library(rworldmap)
# get world map
wmap <- getMap(resolution="high")

firsts_sf <- firsts %>%
  sf::st_as_sf(coords = c("longitude", "latitude"),
               crs = 4326, agr = "constant") 

## Lambert EA: EPSG:3035 ## Lambert CConic: EPSG:3034
wmap_prj <- st_as_sf(wmap) %>% st_transform(crs = "EPSG:3035")
firsts_prj <- firsts_sf %>% st_transform("EPSG:3035")

# bbox_prj <- st_bbox(edge_prj)
bbox_prj <- st_bbox(
  c(xmin = -16, xmax = 21,
    ymin = 8, ymax = 65.2), crs = 4326) %>% 
  st_as_sfc() %>% st_transform("EPSG:3035") %>% st_bbox()

###
map <- ggplot() +
  geom_sf(data = wmap_prj, fill = "grey70", color = "white", size=0.2) +
  # geom_sf(data = arrange(firsts_prj, n_id), 
  #         aes(col = inout), alpha = .5) +
  geom_sf(data = firsts_prj, alpha = .25) +
  coord_sf(xlim = 
             c(bbox_prj[1], bbox_prj[3]), ylim = c(bbox_prj[2], bbox_prj[4]), 
           expand = T) +
  # theme_void() +
  theme(
    plot.background = element_rect(fill = "white"),
    panel.background = element_rect(fill = "white"),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    panel.border = element_rect(colour = "black", fill=NA, size=.5),
    legend.key=element_blank()
  )
map

## SAVE ---------------------------------------------------------------------
ggsave(paste0("figures/trackdeploy_locations.png"),
       plot=map, width=5, height = 6)


##
firsts_bin <- firsts %>% 
  mutate(
    latitude  = round(latitude, 1),
    longitude = round(longitude, 1)
  ) %>% group_by(latitude, longitude) %>% summarise(
    n_id = n_distinct(bird_id),
    n_schemes = n_distinct(scheme)
  )


## 
firsts_bin %>%
  sf::st_as_sf(coords = c("longitude", "latitude"),
               crs = 4326, agr = "constant") %>%
  mapview::mapview()


## Static map: project for prettier map --------------------------------------
library(rworldmap)
# get world map
wmap <- getMap(resolution="high")

firsts_sf <- firsts %>%
  sf::st_as_sf(coords = c("longitude", "latitude"),
               crs = 4326, agr = "constant") 

## Lambert EA: EPSG:3035 ## Lambert CConic: EPSG:3034
wmap_prj <- st_as_sf(wmap) %>% st_transform(crs = "EPSG:3035")
firsts_prj <- firsts_sf %>% st_transform("EPSG:3035")

# bbox_prj <- st_bbox(edge_prj)
bbox_prj <- st_bbox(
  c(xmin = -16, xmax = 21,
    ymin = 8, ymax = 65.2), crs = 4326) %>% 
  st_as_sfc() %>% st_transform("EPSG:3035") %>% st_bbox()

###
map <- ggplot() +
  geom_sf(data = wmap_prj, fill = "grey80", color = "white", size=0.2) +
  # geom_sf(data = arrange(firsts_prj, n_id), 
  #         aes(col = inout), alpha = .5) +
  geom_sf(data = firsts_prj, alpha = .25) +
  coord_sf(xlim = 
             c(bbox_prj[1], bbox_prj[3]), ylim = c(bbox_prj[2], bbox_prj[4]), 
           expand = T) +
  # theme_void() +
  theme(
    plot.background = element_rect(fill = "white"),
    panel.background = element_rect(fill = "white"),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    panel.border = element_rect(colour = "black", fill=NA, size=.5),
    legend.key=element_blank()
  )
map

## SAVE ---------------------------------------------------------------------
ggsave(paste0("figures/cring_locations.png"),
       plot=map, width=5, height = 6)

