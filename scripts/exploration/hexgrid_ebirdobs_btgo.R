## Ebird obs data for Black-tailed godwit

pacman::p_load(auk, mapview, sf, dplyr, terra)
# path to the ebird data file, here a sample included in the package
# get the path to the example data included in the package
# in practice, provide path to ebd, e.g. f_in <- "data/ebd_relFeb-2018.txt
f_in <- "data/ebird/ebd_bktgod_relAug-2022/ebd_bktgod_relAug-2022.txt"
# output text file

f_out <- "data/ebird/filtered/btgo_filtered.txt"

## filter by bounding box
boxll <- st_bbox(c(xmin = -30, xmax = 45, 
                   ymin = 0, ymax = 71.5),
                 crs = 4326)
boxll %>% st_as_sfc() %>% mapview()

## filter input
ebird_data <- f_in %>% 
  auk_ebd() %>%
  # 1. spatial bounding box
  auk_bbox(bbox=boxll) %>% 
  auk_filter(file = f_out, overwrite = T) %>% 
  read_ebd()

## map
ebird_data %>%
  slice(1:38000) %>% 
  select(longitude, latitude) %>%
  st_as_sf(coords = c("longitude", "latitude"),
           crs = 4326, agr = "constant") %>% mapview()

## hexbin the data
library(dggridR)

#Construct a global grid with cells approximately 1000 miles across
# dggs <- dgconstruct(spacing=8, resround='down')
dggs <- dgconstruct(spacing=100, resround='down')

#Get the corresponding grid cells for each bird obs (lat-long pair)
ebird_data$cell <- dgGEO_to_SEQNUM(
  dggs, ebird_data$longitude, ebird_data$latitude)$seqnum

#Get the number of records in each equally-sized cell
binned <- ebird_data %>% group_by(cell) %>% summarise(count=n())

#Get the grid cell boundaries for cells which had obs
grid              <- dgcellstogrid(dggs, ebird_data$cell)
colnames(grid)[1] <- "cell"
cellcenters       <- dgSEQNUM_to_GEO(dggs, grid$cell)

#Update the grid cells' properties to include the number of obs in each cell
grid <- merge(grid, binned, by.x="cell", by.y="cell")

mapview(grid, zcol = "count")
# grid %>% filter(count >1) %>% mapview(zcol = "count")

grid$cc_lat <- cellcenters$lat_deg 

## map it 
ggplot() + 
  geom_sf(data=grid, aes(fill=count), color=NA, alpha=0.4) 

## remove cells in Europe w/ fewer than X obs
grid_f <- grid[which(grid$cc_lat < 35.5 | grid$count > 1),]

mapview(grid_f, zcol="count")

saveRDS(grid_f, "data/analysis/site_hexgrid_ebirdobs_87km.rds")




### ---------------------------------------------------------------------------
## combine hexcells w/ obs that are bordering (NOT WORKING)
sfpolys <- grid

clusterSF <- function(sfpolys, thresh){
  dmat = st_distance(sfpolys)
  hc = hclust(as.dist(dmat>thresh), method="single")
  groups = cutree(hc, h=0.5)
  d = st_sf(
    geom = do.call(c,
                   lapply(1:max(groups), function(g){
                     st_union(sfpolys[groups==g,])
                   })
    )
  )
  d$group = 1:nrow(d)
  d
}

## map empty grid of study area -----------------------------------------------
globgrid <- dggridR::dgearthgrid(dggs)

library(rworldmap)
# get world map
wmap <- getMap(resolution="high")

## Lambert EA: EPSG:3035 ## Lambert CConic: EPSG:3034
wmap_prj <- st_as_sf(wmap) %>% st_transform(crs = "EPSG:3035")
# grid_prj <- globgrid %>% st_transform(crs = "EPSG:3035")
bbox_prj <- st_bbox(
  c(xmin = -16, xmax = 21,
    ymin = 8, ymax = 65.2), crs = 4326) %>% 
  st_as_sfc() %>% st_transform("EPSG:3035") %>% st_bbox()

###
map <- ggplot() +
  geom_sf(data = wmap_prj, fill = "grey60", color = "white", size=0.3) +
  geom_sf(data = globgrid, 
          aes(), fill = NA, col = "black", size=0.1) + #, alpha = 0.65
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

## SAVE map

ggsave(paste0("figures/empty_hex87km.png"), plot=map, width=5, height = 6)


## Map ebird areas w/ BTGO sightings -------------------------------
# cit sci grid
grid_cs <- readRDS("data/analysis/site_hexgrid_ebirdobs_87km.rds")
# grid_cs_prj <- st_as_sf(grid_cs) %>% st_transform(crs = "EPSG:3035")


###
map2 <- ggplot() +
  geom_sf(data = wmap_prj, fill = "grey70", color = NA) +
  geom_sf(data = grid_cs,
          aes(fill = count), color=NA) +
  geom_sf(data = wmap_prj, fill = NA, color = "white", size=0.2) +
  coord_sf(xlim =
             c(bbox_prj[1], bbox_prj[3]), ylim = c(bbox_prj[2], bbox_prj[4]),
           expand = T) +
  theme(
    plot.background = element_rect(fill = "white"),
    panel.background = element_rect(fill = "white"),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    panel.border = element_rect(colour = "black", fill=NA, size=.5),
    legend.key=element_blank()
  ) +
  labs(fill='Checklists') 

ggsave(paste0("figures/hex87km_citscisitesX.png"), plot=map2, width=5, height = 6)
