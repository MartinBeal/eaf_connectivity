## Create a spatial, undirected network from ring recapture and resighting data
# to sites defined using a polygon layer (IBAs, RAMSAR etc)

pacman::p_load(dplyr, igraph, stringr, tictoc, tidygraph, sfnetworks, ggplot2, 
               sf, mapview, magrittr, lubridate, netrankr)


## Run through each data type -------------------------------------------------
# datatype <- "metal"
datatype <- "color"
# datatype <- "trax"

if(datatype == "color"){
  ## color ringed bird captures and resightings overlaid on polygon layer
  alldat <- readRDS("data/analysis/ringing/cring_merge_no7dayreobs_ibas.rds")
  alldat %<>% rename(timestamp = date)
} else if (datatype == "metal"){
  ## metal ring captures, recaptures, recoveries
  alldat <- readRDS("data/analysis/ringing/euring_metal_ibas.rds")
  alldat %<>% rename(timestamp = date)
} else if (datatype == "trax"){
  ## tracking locations 
  # alldat <- readRDS("data/analysis/tracking/PTT_GPS_mconn_12h_ibas.rds")
  alldat <- readRDS("data/analysis/tracking/PTT_GPS_mconn_12h_no0_ibas.rds")
  # alldat %<>% rename(timestamp = date)
}

# fxn for splitting string into columns
source("C:/Users/Martim Bill/Documents/R/source_scripts/str2col.R")
# fxn for calc. displacement distance from first tracking locations
source("C:/Users/Martim Bill/Documents/R/source_scripts/displ_dist.R")

# alldat %<>% rename(combo = metal)
# alldat %<>% filter(obssource == "euring") # which data sources
# alldat %<>% filter(obssource != "euring")

## choose data subset to run --------------------------------------------------

## which network to create
# season <- "all"
# season <- "spring"
season <- "fall"

### Separate networks for fall and spring migration
## Spring: January 1 - June 30th, Fall: June 24th - January 31st

if(season == "all"){ ## year-round
  netdat <- alldat ## all individuals
} else if(season == "spring"){
  netdat <- subset(alldat, month(alldat$timestamp) %in% c(1:6))
} else if(season == "fall"){
  doy <- lubridate::yday(alldat$timestamp) # June 23/24: 175
  netdat <- subset(
    alldat, 
    month(alldat$timestamp) %in% c(1, 7:12) | doy %in% c(175:181)
  )
}

## Tracking data: rmv ids w/ only local displacement in season (i.e. no mig.)
if(datatype == "trax"){
  netdat <- rename(netdat, id = bird_id)
  ## add column of diplacement distance from first location 
  netdat <- displ_dist(netdat)
  
  ## Summarise max/avg distance from first points -------------------------------
  
  displ_id <- netdat %>% group_by(id) %>% 
    summarise(
      mn_displ = mean(disp_km),
      sd_displ = sd(disp_km),
      md_displ = median(disp_km),
      mx_displ = max(disp_km)
    )
  
  ## remove individuals w/ only local displacement (i.e. no migration)
  
  wmig <- filter(displ_id, mx_displ >= 200)
  
  n_distinct(netdat$id)
  netdat %<>% filter(id %in% wmig$id)
  n_distinct(netdat$id) 
  
  netdat <- rename(netdat, bird_id = id)
  
}

## show data from a certain place
# netdat <- subset(netdat, bird_id %in% unique(alldat$bird_id)[1:100])
# netdat <- subset(netdat, scheme_country == "Iceland")

## top sites visited
# xx <- netdat %>% group_by(SitRecID, IntName, country) %>% 
#   summarise(
#     n_birds = n_distinct(bird_id),
#     n_obs   = n()
#   ) %>% arrange(desc(n_birds)) %>% ungroup() %>% 
#   filter(country != "Iceland") %>% 
#   slice(1:11)

## 
netdat %<>% rename(site_poly = IntName) ## IBAs only


###---------------------------------------------------------------------------
### network ------------------------------------------------------------------

## (optionally) remove locations falling outside any site polygon
nones  <- filter(netdat, site_poly == "none")

## filter to only obs w/in sites (and for trax, only site visits >= 2 days)
if(datatype == "color"){
  saveRDS(
    nones, 
    paste0("data/analysis/ringing/color_outside_", season, "_ibas10km.rds"))
  netdat <- filter(netdat, site_poly != "none")
} else if (datatype == "metal"){
  saveRDS(
    nones, 
    paste0("data/analysis/ringing/metal_outside_", season, "_ibas10km.rds"))
  netdat <- filter(netdat, site_poly != "none")
} else if (datatype == "trax"){
  saveRDS(
    nones, 
    paste0("data/analysis/tracking/outside_trax_", season, "_ibas10km.rds"))
  netdat <- filter(netdat, site_poly != "none" & n_day >= 2)
}

## again remove birds w/ only one sighting (to avoid unconnected nodes) -------
##**WILL NEED TO ADD BACK IN TO MAKE COMPARISONS W/ OUTSIDE POINTS ##
xz <- netdat %>% group_by(bird_id) %>% 
  summarise(nobs = n()) %>% filter(nobs == 1)
netdat <- subset(netdat, !bird_id %in% xz$bird_id)

## remove birds only seen (multiple times) at same site ----------------------- 
##**ALSO HERE WILL NEED TO ADD BACK IN TO MAKE COMPARISONS W/ OUTSIDE POINTS ##
nsites <- netdat %>% group_by(bird_id) %>% 
  summarise(nsites = n_distinct(site_poly))

## % of birds w/ relocs at one site only
sum(nsites$nsites == 1) / n_distinct(netdat$bird_id) * 100

xz2 <- filter(nsites, nsites == 1)
netdat <- subset(netdat, !bird_id %in% xz2$bird_id)

## (local) numeric code of nodes/sites ------------------------------------------------
# netdat$loc_num <- as.numeric(as.factor(netdat$site_poly)) # name
netdat$loc_num <- as.numeric(as.factor(netdat$SitRecID))  # (absolute) numeric

# create reference table of sites w/ obs --------------------------------------
site_cent <- readRDS( ## site centroids
  "data/geodata/ibas/Africa_Europe_IBA/Africa_Europe_IBA_centroids.shp") 
site_cent %<>% 
  mutate(SitRecID = as.character(SitRecID)) %>% 
  dplyr::select(SitRecID, geometry)

site_summ <- netdat %>% 
  group_by(loc_num, site_poly, SitRecID) %>% summarise() %>% 
  left_join(site_cent, by = c("SitRecID")) %>% st_as_sf()

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
n_id_total <- n_distinct(noself$bird_id)

edgelist <- noself %>% group_by(sitecomb) %>% 
  summarise(
    n_id = n(),
    prop_id = n_id / n_id_total
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
  dplyr::select(from, to, n_id, prop_id) %>% 
  as_tibble()


## Vertex list ---------------------------------------------------------------
if(datatype %in% c("color", "metal")){
  n_obstype <- netdat %>% 
    group_by(loc_num, SitRecID, obstype) %>% 
    summarise(
      n_obs = n()
    ) %>% tidyr::pivot_wider(
      # id_cols = "obstype",
      names_from = "obstype",
      names_prefix = "n_",
      values_from = "n_obs"
    )
} else if (datatype == "trax"){
  n_obstype <- netdat %>% 
    group_by(loc_num, SitRecID, device) %>% 
    summarise(
      n_obs = n()
    ) %>% tidyr::pivot_wider(
      names_from = "device",
      names_prefix = "n_",
      values_from = "n_obs"
    )
}

n_id_total <- n_distinct(netdat$bird_id)

nodelist <- netdat %>% group_by(loc_num, SitRecID,) %>% 
  summarise(
    n_id    = n_distinct(bird_id),
    prop_id = n_id / n_id_total,
    n_obs   = n()
  ) %>% right_join(n_obstype)

nodelist <- nodelist %>% 
  left_join(site_summ) %>% 
  sf::st_as_sf()

## Convert to sfnetwork -------------------------------------------------------

# netsf <- sfnetwork(nodelist, edgelist, node_key = "dist", directed = F)
netsf <- sfnetwork(nodelist, edgelist, node_key = "dist", directed = F, 
                   edges_as_lines = TRUE)


## filter out sites (and links) w/ few observations
# netsf <- netsf %>%
#   activate("nodes") %>%
#   filter(n_id > 1)


### Calculate network metrics ------------------------------------------------

# netsf <- readRDS("data/analysis/networks/metal_all_iba10km_poly.rds") # metal
# netsf <- readRDS("data/analysis/networks/color_all_iba10km_poly.rds") # color
# netsf <- readRDS("data/analysis/networks/trax_all_iba10km_poly.rds")  # trax

# Degree = the number of adjacent edges for a node
# Betweenness = the number of shortest paths going through a node

## global metrics 
netsize <- nrow(st_as_sf(netsf, "nodes")) # network size (n nodes)
nedges  <- nrow(st_as_sf(netsf, "edges"))

## ego metrics (node/edge-level)
netsf %<>%
  activate(nodes) %>%
  mutate(
    degree      = centrality_degree(),
    degree_norm = degree / n_distinct(loc_num), # normalized 0-1 (prop of sites)
    degree_rank = dense_rank(desc(degree)),
    between     = centrality_betweenness(directed = FALSE), # node betweenness
    between_norm = centrality_betweenness(directed = FALSE, normalized = T), # normalized
    btwn_rank   = dense_rank(desc(between))
  )

## interactive map it
nodesf <- netsf %>% activate("nodes") %>% sf::st_as_sf()
edgesf <- netsf %>% activate("edges") %>% sf::st_as_sf()
#   mapview::mapview(nodesf, zcol="n_id")

## just nodes
# mapview::mapview(nodesf, zcol="n_id")
# mapview::mapview(nodesf, zcol="between")
# mapview::mapview(nodesf, zcol="degree")
# mapview::mapview(nodesf, zcol="between_norm")
mapview::mapview(nodesf, zcol="degree_rank")
# mapview::mapview(nodesf, zcol="btwn_rank")

# mapview::mapview(nodesf, zcol="btwn_rank") 
# mapview::mapview(nodesf, zcol="btwn_rank") +
# (filter(edgesf, from %in% 241 | to %in% 241) %>% 
#   mapview::mapview())

## SAVE ##

saveRDS(netsf, paste0("data/analysis/networks/", datatype,"_", season, "_iba10km_poly.rds"))

