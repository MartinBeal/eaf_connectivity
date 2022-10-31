## Identify 'sites' used by tracked birds

pacman::p_load(dggridR, magrittr, sfnetworks, ggplot2, data.table)

td <- readRDS("data/analysis/tracking/PTT_GPS_mconn_12h.rds")

# fxn for splitting string into columns
source("C:/Users/Martim Bill/Documents/R/source_scripts/str2col.R")

## choose data subset to run --------------------------------------------------

## which network to create
season <- "all"
# season <- "spring"
# season <- "fall"

# netdat <- subset(alldat, bird_id %in% unique(alldat$bird_id)[1:100])

if(season == "all"){ ## year-round
  nettd <- td ## all individuals
} else { ## separate networks for fall and spring migration 
  td$season <- ifelse(month(td$timestamp) < 7, "spring", "fall")
  nettd <- td[which(td$season == season), ]
}

## hexbin obs per individual -- end up w/ only points where bird is stopped
# hexcells where more than one day spent

#Construct a global grid with cells approximately 100 km across
# dggs <- dgconstruct(spacing=8, resround='down')
dggs <- dgconstruct(spacing=100, resround='down')

## create hexgrid from td obs
nettd$cell <- dgGEO_to_SEQNUM(
  dggs, nettd$longitude, nettd$latitude)$seqnum

## Get number of days spent in each cell per individual
nettd$year <- year(nettd$timestamp)
# binned_td <- nettd %>% group_by(id, year, cell) %>% 
#   summarise(N = n())
binned_td <- as.data.table(nettd)[, .N, by=list(id, year, cell)]

# ## remove cells w/ only one point (<=12 h)
binned_td <- filter(binned_td, N > 1)

## add up points/days across birds
bin_agg <- binned_td %>% group_by(cell) %>% 
  summarise(
    n_locs = sum(N),
    n_birds = n_distinct(id))


#Get the grid cell boundaries for cells which had obs
grid_td           <- dgcellstogrid(dggs, bin_agg$cell)
colnames(grid_td)[1] <- "cell"

## Get cell centers (to be nodes)
ccenters_td       <- dgSEQNUM_to_GEO(dggs, bin_agg$cell)

#Update the grid cells' properties to include the number of obs in each cell
grid_td <- merge(grid_td, bin_agg, by.x="cell", by.y="cell")
 
# grid_td %>% # map
#   sf::st_as_sf(
#     coords = c("lon", "lat"),
#     crs = 4326, agr = "constant") %>% mapview(zcol="n_locs")


## add cell center coordinates 
grid_td <- grid_td %>% st_drop_geometry() %>% 
  bind_cols(latitude=ccenters_td$lat_deg, longitude=ccenters_td$lon_deg)


###---------------------------------------------------------------------------
### network ------------------------------------------------------------------

## filter tracking data to cells w/ more than X amount of relocs
nettd_f <- filter(nettd, cell %in% unique(binned_td$cell))

## numeric coding of locations
nettd_f$loc_num <- as.numeric(as.factor(nettd_f$cell))

nettd_f <- nettd_f %>% 
  dplyr::select(-latitude, -longitude) %>% 
  left_join(grid_td, by="cell")

# create reference table of sites w/ obs
site_summ <- nettd_f %>% group_by(loc_num, cell) %>% 
  summarise(latitude = first(latitude), longitude = first(longitude))


## Edge list ------------------------------------------------------------------

## produce all combinations of sites visited by each individual
nettd_list <- split(nettd_f, nettd_f$id)

tictoc::tic()
nettd_list <- lapply(
  seq_along(nettd_list), 
  function(x){
    # print(x)
    one <- nettd_list[[x]]
    xx <- as.data.frame(
      RcppAlgos::comboGrid(
        one$loc_num,
        one$loc_num, repetition = TRUE)
    )
    xx$id <- one$id[1]
    return(xx)
  })
tictoc::toc()

full <- data.table::rbindlist(nettd_list)

## remove self connections
noself <- full[-which(full$Var1 == full$Var2), ]

## combine sites into single variable for summarizing
noself$sitecomb <- paste(noself$Var1, noself$Var2)

## retain only unique site combos, summ how many individuals for connex
edgelist <- noself %>% group_by(sitecomb) %>% 
  summarise(
    n_id = n()
  )

## re-split site combos into separate columns
edgelist <- cbind(
  edgelist,
  str2col(edgelist$sitecomb,
          pattern = " ", 
          cols = 1:2, colnames = c("from", "to"))
) %>% 
  mutate(
    from = as.integer(from), to = as.integer(to)
  ) %>% 
  dplyr::select(from, to, n_id) %>% 
  as_tibble()



## Vertex list ---------------------------------------------------------------
# n_obstype <- nettd_f %>% 
  # group_by(loc_num) %>% 
  # summarise(
    # n_obs = n()
  # )# %>% tidyr::pivot_wider(
    # id_cols = "obstype",
    #names_from = "obstype",
    #names_prefix = "n_",
    #values_from = "n_obs"
  #)

nodelist <- nettd_f %>% group_by(loc_num) %>% 
  summarise(
    n_id = n_distinct(id),
    n_obs  = n()
  ) # %>% right_join(n_obstype)

nodelist <- nodelist %>% 
  left_join(site_summ) %>% 
  sf::st_as_sf(coords = c("longitude", "latitude"), 
               crs = 4326, agr = "constant")

## Convert to sfnetwork -------------------------------------------------------

# netsf <- sfnetwork(nodelist, edgelist, node_key = "dist", directed = F)
netsf <- sfnetwork(nodelist, edgelist, node_key = "dist", directed = F, 
                   edges_as_lines = TRUE)

## if it gives an error try running via igraph
# graph_from_data_frame(
#   d = edgelist, 
#   vertices = nodelist,
#   directed = F)

## filter out sites (and links) w/ few observations
# netsf <- netsf %>%
#   activate("nodes") %>%
#   filter(n_id > 1)

## interactive map it
nodesf <- netsf %>% activate("nodes") %>% sf::st_as_sf()
# netsf %>% activate("edges") %>% sf::st_as_sf() %>% mapview::mapview() +
#   mapview::mapview(nodesf, zcol="n_id")
## just nodes
# mapview::mapview(nodesf, zcol="n_id")


## Static map: project for prettier map --------------------------------------
## get country centroids 
library(rworldmap)
# get world map
wmap <- getMap(resolution="high")

## Lambert EA: EPSG:3035 ## Lambert CConic: EPSG:3034
wmap_prj <- st_as_sf(wmap) %>% st_transform(crs = "EPSG:3035")
edge_prj <- st_as_sf(netsf, "edges") %>% st_transform(crs = "EPSG:3035")
node_prj <- st_as_sf(netsf, "nodes") %>% st_transform("EPSG:3035")
bbox_prj <- st_bbox(
  c(xmin = -16, xmax = 21,
    ymin = 8, ymax = 65.2), crs = 4326) %>% 
  st_as_sfc() %>% st_transform("EPSG:3035") %>% st_bbox()

###
map <- ggplot() +
  geom_sf(data = wmap_prj, fill = "grey80", color = NA) +
  geom_sf(data = edge_prj,
          aes(size = n_id), col = "black") + #, alpha = 0.65
  geom_sf(data = wmap_prj, fill = NA, color = "white", size=0.2) +
  geom_sf(data = node_prj, 
          aes(col = n_id)) +
  coord_sf(xlim = 
             c(bbox_prj[1], bbox_prj[3]), ylim = c(bbox_prj[2], bbox_prj[4]), 
           expand = T) +
  scale_size(range = c(0.01, 2)) +
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

ggsave(paste0("figures/", season, "_traxX.png"), plot=map, width=5, height = 6)
# ggsave(paste0("figures/", season, "_trax_nodes.png"), plot=map, width=5, height = 6)


## ----------------------------------------------------------------------------
## combine map w/ ebird areas w/ BTGO sightings -------------------------------
# cit sci grid
grid_cs <- readRDS("data/analysis/site_hexgrid_ebirdobs_87km.rds")
# grid_cs_prj <- st_as_sf(grid_cs) %>% st_transform(crs = "EPSG:3035")

###
map2 <- ggplot() +
  geom_sf(data = wmap_prj, fill = "grey80", color = NA) +
  geom_sf(data = grid_cs,
          aes(fill = count), color=NA) +
  # geom_sf(data = edge_prj,
  #         aes(size = n_id), color = "black") + #, alpha = 0.65
  geom_sf(data = wmap_prj, fill = NA, color = "white", size=0.2) +
  geom_sf(data = node_prj,
          color="red", size=.5) +
  coord_sf(xlim =
             c(bbox_prj[1], bbox_prj[3]), ylim = c(bbox_prj[2], bbox_prj[4]),
           expand = T) +
  scale_size(range = c(0.01, 2)) +
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

ggsave(paste0("figures/", season, "_trax_hex87km_citscisitesX.png"), plot=map2, width=5, height = 6)
