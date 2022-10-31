## Combine in-site and out-site nodes

outsite <- nodesf
# insite  <- nodesf

outsite$loc <- paste0("o_", outsite$loc_num)
insite$loc  <- paste0("i_", insite$loc_num)

outsite$inout <- "out"
insite$inout  <- "in"

allsites <- bind_rows(
  outsite[,c("loc", "inout", "loc_num", "n_id", "n_obs", "geometry")],
  insite[,c("loc", "inout", "loc_num", "n_id", "n_obs", "geometry")])

allsites %>% mapview::mapview(zcol="inout")


## Static map: project for prettier map --------------------------------------
library(rworldmap)
# get world map
wmap <- getMap(resolution="high")

## Lambert EA: EPSG:3035 ## Lambert CConic: EPSG:3034
wmap_prj <- st_as_sf(wmap) %>% st_transform(crs = "EPSG:3035")
node_prj <- allsites %>% st_transform("EPSG:3035")

# bbox_prj <- st_bbox(edge_prj)
bbox_prj <- st_bbox(
  c(xmin = -16, xmax = 21,
    ymin = 8, ymax = 65.2), crs = 4326) %>% 
  st_as_sfc() %>% st_transform("EPSG:3035") %>% st_bbox()

###
map <- ggplot() +
  geom_sf(data = wmap_prj, fill = "grey80", color = NA) +
  # geom_sf(data = edge_prj,
  #         aes(size = n_id), col = "black") + #, alpha = 0.65
  geom_sf(data = wmap_prj, fill = NA, color = "white", size=0.2) +
  geom_sf(data = arrange(node_prj, n_id), 
          aes(col = inout), alpha = .5) +
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
# map

## SAVE ---------------------------------------------------------------------
ggsave(paste0("figures/", datatype,"_", season, "_inandout_iba10kmX.png.png"),
       plot=map, width=5, height = 6)
