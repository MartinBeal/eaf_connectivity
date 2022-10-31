## Create a spatial, undirected network from ring recapture and resighting data

pacman::p_load(dplyr, igraph, stringr, tictoc, tidygraph, sfnetworks, ggplot2, 
               sf, mapview, magrittr, lubridate)

# alldat <- readRDS("data/analysis/test_comb_euring_cring_districts.rds") # county level
## ice, wash, tagus, dutch, and EURING data:
alldat <- readRDS("data/analysis/ringing/test_comb_euring_cring_no7dayreobs.rds")


# fxn for splitting string into columns
source("C:/Users/Martim Bill/Documents/R/source_scripts/str2col.R")

# alldat %<>% rename(combo = metal)

# alldat %<>% filter(obssource == "euring") # which data sources
# alldat %<>% filter(obssource != "euring")

## choose data subset to run --------------------------------------------------

# level <- "site"
# level <- "district"
# level <- "country"

## which network to create
season <- "all"
# season <- "spring"
# season <- "fall"

# netdat <- subset(alldat, bird_id %in% unique(alldat$bird_id)[1:100])

if(season == "all"){ ## year-round
  netdat <- alldat ## all individuals
} else { ## separate networks for fall and spring migration 
  alldat$season <- ifelse(month(alldat$date) < 7, "spring", "fall")
  netdat <- alldat[which(alldat$season == season), ]
}

## only Euring data (inspection)
# netdat <- subset(netdat, obssource == "euring")

### Spatial scale -------------------------------------
# netdat <- alldat
## for starters, round lat/longs to X decimal places to aggregate obs into sites

# netdat$latitude  <- round(netdat$latitude, 2)
# netdat$longitude <- round(netdat$longitude, 2)

## numeric coding of locations
netdat$loc_num <- as.numeric(as.factor(paste(netdat$latitude, netdat$longitude)))

# create reference table of sites w/ obs
site_summ <- netdat %>% group_by(loc_num) %>% 
  summarise(latitude = first(latitude), longitude = first(longitude))

## Quickly get country centroids 
library(rgeos)
library(rworldmap)
# get world map
wmap <- getMap(resolution="high")


## Edge list ------------------------------------------------------------------

## produce all combinations of sites visited by each individual
netdat_list <- split(netdat, netdat$bird_id)

tic()
netdat_list <- lapply(
  seq_along(netdat_list), 
  function(x){
    # print(x)
    one <- netdat_list[[x]]
    xx <- as.data.frame(
      RcppAlgos::comboGrid(
        one$loc_num,
        one$loc_num, repetition = TRUE)
    )
    xx$bird_id <- one$bird_id[1]
    return(xx)
  })
toc()

full <- data.table::rbindlist(netdat_list)

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
n_obstype <- netdat %>% 
  group_by(loc_num, obstype) %>% 
  summarise(
    n_obs = n()
  ) %>% tidyr::pivot_wider(
    # id_cols = "obstype",
    names_from = "obstype",
    names_prefix = "n_",
    values_from = "n_obs"
  )

nodelist <- netdat %>% group_by(loc_num) %>% 
  summarise(
    n_id = n_distinct(bird_id),
    n_obs  = n()
  ) %>% right_join(n_obstype)

nodelist <- nodelist %>% 
  left_join(site_summ) %>% 
  sf::st_as_sf(coords = c("longitude", "latitude"), 
               crs = 4326, agr = "constant")

## Convert to sfnetwork

# netsf <- sfnetwork(nodelist, edgelist, node_key = "dist", directed = F)
netsf <- sfnetwork(nodelist, edgelist, node_key = "dist", directed = F, 
                   edges_as_lines = TRUE)

## if it gives an error try running via igraph
# graph_from_data_frame(
#   d = edgelist, 
#   vertices = nodelist,
#   directed = F)

## filter out sites (and links) w/ few observations
netsf <- netsf %>%
  activate("nodes") %>%
  filter(n_id > 1)

## interactive map it
nodesf <- netsf %>% activate("nodes") %>% sf::st_as_sf()
# netsf %>% activate("edges") %>% sf::st_as_sf() %>% mapview::mapview() +
#   mapview::mapview(nodesf, zcol="n_id")
## just nodes
mapview::mapview(nodesf, zcol="n_id")


## Static map: project for prettier map --------------------------------------

## Lambert EA: EPSG:3035 ## Lambert CConic: EPSG:3034
wmap_prj <- st_as_sf(wmap) %>% st_transform(crs = "EPSG:3035")
edge_prj <- st_as_sf(netsf, "edges") %>% st_transform(crs = "EPSG:3035")
node_prj <- st_as_sf(netsf, "nodes") %>% st_transform("EPSG:3035")
bbox_prj <- st_bbox(edge_prj)

###
map <- ggplot() +
  geom_sf(data = wmap_prj, fill = "grey80", color = NA) +
  geom_sf(data = edge_prj, 
          aes(size = n_id), col = "black") + #, alpha = 0.65
  geom_sf(data = wmap_prj, fill = NA, color = "white") +
  geom_sf(data = node_prj, 
          aes(col = n_id), size=0.2) +
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

ggsave(paste0("figures/", season, "_latlondec1_tag_wash_ice.png"), plot=map, width=5, height = 6)


### Calculate network metrics ------------------------------------------------

library(netrankr)

netsf %<>%
  activate(nodes) %>%
  mutate(
    between = centrality_betweenness() # node betweenness
  )

ggplot() +
  geom_sf(data = st_as_sf(wmap), fill = "grey90", color = "grey70") +
  # geom_sf(data = st_as_sf(netsf, "edges"), col = "grey30") +
  geom_sf(data = st_as_sf(netsf, "nodes"), aes(color = between, size = between)) +
  coord_sf(xlim = c(bbox[1], bbox[3]), ylim = c(bbox[2], bbox[4]), expand = T) +
  theme_void() +
  theme(
    panel.background = element_rect(fill = "white") 
  )
