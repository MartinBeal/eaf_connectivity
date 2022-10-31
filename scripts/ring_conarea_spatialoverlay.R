## Overlap relocations on polygon layer (Ramsar/IBA) 

pacman::p_load(magrittr, dplyr, ggplot2, sf, rworldmap)

rams <- raster::shapefile("data/geodata/ramsar/EAF_ramsar.shp")
iba <- raster::shapefile("data/geodata/ibas/EAF_btgo_IBA/eaf_btgo_iba.shp")


## map ibas
# get world map
wmap <- getMap(resolution="high")

rams_sf_prj <- st_as_sf(rams) %>% st_transform(crs = "EPSG:3035")
iba_sf_prj <- st_as_sf(iba) %>% st_transform(crs = "EPSG:3035")
wmap_prj <- st_as_sf(wmap) %>% st_transform(crs = "EPSG:3035")
bbox_prj <- st_bbox( ## manual limits from network data (change as data changes)
  c(xmin = 1271137.2, xmax = 5214780.5, ymax = 5133771.9, ymin = -720268.5),
  crs = "EPSG:3035")

###
map <- ggplot() +
  geom_sf(data = wmap_prj, aes(), fill = "white") +
  geom_sf(data = iba_sf_prj, aes(), fill = "red", color=NA) +
  geom_sf(data = rams_sf_prj, aes(), fill = "blue", color=NA) +
  coord_sf(xlim = 
             c(bbox_prj[1], bbox_prj[3]), ylim = c(bbox_prj[2], bbox_prj[4]), 
           expand = T) +
  scale_fill_continuous(low="#fee0d2", high="#de2d26", 
                        guide="colorbar", na.value="white") +
  # theme_void() +
  theme(
    plot.background = element_rect(fill = "white"),
    # panel.background = element_rect(fill = "slategray1"),
    panel.background = element_rect(fill = "grey85"),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    panel.border = element_rect(colour = "black", fill=NA, size=.5),
    legend.key=element_blank()
  )
map

ggsave("figures/iba_ramsar.png", height = 7, width=5)

## Identify closest polygon for points not overlapping any

# Iterate over all points using sapply
sapply(1:nrow(points_sf), function(x) min(st_distance(roads_sf, points_sf[x, ])))