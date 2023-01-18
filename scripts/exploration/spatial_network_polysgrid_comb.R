### Create a spatial, undirected network from ring recapture and resighting data
## to sites defined using a polygon layer (IBAs, RAMSAR etc)
# VERSION that combines in and out sites into network

pacman::p_load(dplyr, igraph, stringr, tictoc, tidygraph, sfnetworks, ggplot2, 
               sf, mapview, magrittr, lubridate, netrankr)


## Run through each data type -------------------------------------------------
# datatype <- "metal"
# datatype <- "color"
datatype <- "trax"

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

##@ Create hexgrid cell 'sites' for points falling outside known site layer ---
## separate locations falling outside any site polygon
outdat  <- filter(netdat, site_poly == "none")
indat  <- filter(netdat, site_poly != "none")

## filter to only obs w/in sites (and for trax, only site visits >= 2 days)
# if(datatype == "color"){
#   saveRDS(
#     outdat, 
#     paste0("data/analysis/ringing/color_outside_", season, "_ibas10km.rds"))
# } else if (datatype == "metal"){
#   saveRDS(
#     outdat, 
#     paste0("data/analysis/ringing/metal_outside_", season, "_ibas10km.rds"))
# } else if (datatype == "trax"){
#   saveRDS(
#     outdat, 
#     paste0("data/analysis/tracking/outside_trax_", season, "_ibas10km.rds"))
# }


## Hexgrid --------------------------------------------------------------------

## Construct global grid with cells approximately 1000 miles across
# dggs <- dgconstruct(spacing=8, resround='down')
dggs <- dgconstruct(spacing=10, resround='down')

## extract cell # for each relocation
outdat$cell <- dgGEO_to_SEQNUM(
  dggs, outdat$longitude, outdat$latitude)$seqnum

outdat$SitRecID <- paste0("cell_", as.character(outdat$cell)) # use cell # as absolute ID

## re-merge in and out data
netdat <- indat %>% mutate(cell = NA) %>% bind_rows(outdat) %>%
  arrange(bird_id, timestamp)


## Calculate number of consecutive days spent at a site -----------------------
netdat <- as.data.table(netdat)
if(datatype == "trax"){
  netdat[, c("samesite","n_day"):=NULL] # remove cols
}

## id consecutive obs at a site/cell, per bird (data.table way)
netdat[, samesite := data.table::rleid(SitRecID), by = bird_id]
netdat$samesite <- ifelse(is.na(netdat$SitRecID), NA, netdat$samesite) # NA where site NA

## number of consecutive days at a site
netdat <- netdat %>%
  group_by(bird_id, samesite) %>%
  summarise(n_day = n_distinct(yday(timestamp))) %>%
  left_join(netdat)

## for tracking data, remove data from sites visited for < 48h
if (datatype == "trax"){
  netdat <- filter(netdat, n_day >= 2)
}

## check out some data
# netdat %>% filter(bird_id == "B.tenskar") %>% 
# netdat[5000:10000,] %>%
#   sf::st_as_sf(
#     coords = c("longitude", "latitude"),
#     crs = 4326, agr = "constant") %>% mapview(zcol="SitRecID")

## Remove any birds w/ only one sighting (to avoid unconnected nodes) ---------
xz <- netdat %>% group_by(bird_id) %>% 
  summarise(nobs = n()) %>% filter(nobs == 1)
netdat <- subset(netdat, !bird_id %in% xz$bird_id)

## remove birds only seen (multiple times) at same site ----------------------- 
nsites <- netdat %>% group_by(bird_id) %>% 
  summarise(nsites = n_distinct(SitRecID))

## % of birds w/ relocs at one site only
sum(nsites$nsites == 1) / n_distinct(netdat$bird_id) * 100

xz2 <- filter(nsites, nsites == 1)
netdat <- subset(netdat, !bird_id %in% xz2$bird_id)

### DOUBLE CHECK THAT THIS STEP IS NEEDED HERE OR IF ONLY LATER WORKS
## create (relative) numeric code for nodes/sites -----------------------------
# netdat$loc_num <- as.numeric(as.factor(netdat$site_poly)) # name
# netdat$loc_num <- as.numeric(as.factor(netdat$SitRecID))  # (absolute) numeric

## split data again by in and out site to add cell center coords --------------
outdat <- netdat %>% 
  filter(site_poly == "none")
indat <- netdat %>% 
  filter(site_poly != "none")

## one row per cell
binned_re <- outdat %>% group_by(SitRecID, cell) %>% 
  summarise() # %>% summarise(
# count=n(), 
# n_birds = n_distinct(bird_id))

## Get the grid cell boundaries for cells which had obs
grid_re           <- dgcellstogrid(dggs, binned_re$cell)
colnames(grid_re)[1] <- "cell"

ccenters_re       <- dgSEQNUM_to_GEO(dggs, grid_re$cell)

## interactive map of cell centers
# cbind.data.frame(lat=ccenters_re$lat_deg, lon=ccenters_re$lon_deg) %>%
#   sf::st_as_sf(
#     coords = c("lon", "lat"),
#     crs = 4326, agr = "constant") %>% mapview()

### Update the grid cells' properties to include the number of obs in each cell
# grid_re <- merge(grid_re, binned_re, by.x="cell", by.y="cell")
# mapview(grid_re, zcol = "count")
# mapview(grid_re, zcol = "n_birds")


## ----------------------------------------------------------------------------
## Group neighbouring hexcells into 'sites' -----------------------------------

## solution found here:
## https://gis.stackexchange.com/questions/323038/dissolve-only-overlapping-polygons-in-r-using-sf

## union intersecting polygons into combined polygon
parts <- st_cast(st_union(grid_re),"POLYGON")
plot(parts)

## identify which 'part' each cell polygon intersects
cellgrp <- unlist(st_intersects(grid_re, parts))

## union cell polygons by their cellgroup membership + paste cell numbers
grid_re2 <- cbind(grid_re, cellgrp) %>%
  group_by(cellgrp) %>%
  summarize(
    cell = paste(cell, collapse = ", ")
  )

mapview(grid_re2)

## centroid of cell group polygon
ccenters_re <- st_centroid(grid_re2) %>% 
  dplyr::mutate(lon_deg = sf::st_coordinates(.)[,1],
                lat_deg = sf::st_coordinates(.)[,2])

## interactive map of cell centers
# cbind.data.frame(lat=ccenters_re$lat_deg, lon=ccenters_re$lon_deg) %>%
#   sf::st_as_sf(
#     coords = c("lon", "lat"),
#     crs = 4326, agr = "constant") %>% mapview()

### Update the grid cells' properties to include the number of obs in each cell
# grid_re2 <- merge(grid_re2, binned_re, by.x="cell", by.y="cell")
# mapview(grid_re, zcol = "count")
# mapview(grid_re, zcol = "n_birds")

## add cell center coordinates
grid_re2 <- grid_re2 %>% 
  bind_cols(latitude=ccenters_re$lat_deg, longitude=ccenters_re$lon_deg)


### Overlay relocations on cell group polygons --------------------------------

outdat %<>% sf::st_as_sf(
  coords = c("longitude", "latitude"),
  crs = 4326, agr = "constant")

ov <- sapply(
  st_intersects(outdat, grid_re2),
  function(x){
    if(length(x) == 0){x <- 'none'}
    return(x[1])
  })

grid_re2$rowid <- 1:nrow(grid_re2)

## combine site info with overlap result
pntscellgrps <- left_join(data.frame(rowid = ov), st_drop_geometry(grid_re2))

## combine overlap result back into bird locations w/ cellgroup info
outdat <- bind_cols(outdat, pntscellgrps[,c("cellgrp", "cell")]) %>% 
  mutate(
    cellgrp = paste0("cellgrp_", cellgrp),
    SitRecID = cellgrp ## make the cell group the main site ID
  )

## create reference table of sites w/ obs --------------------------------------
## Out site centroids:
out_cent <- ccenters_re %>% 
  mutate(
    cellgrp = paste0("cellgrp_", cellgrp),
    SitRecID = cellgrp ## make the cell group the main site ID
  ) %>% group_by(SitRecID, cellgrp) %>% 
  st_as_sf(coords = c("longitude", "latitude"), 
           crs = 4326, agr = "constant") %>% 
  summarise() %>% 
  dplyr::select(SitRecID, geometry)


## centroids of site polygons
site_cent <- readRDS( 
  "data/geodata/ibas/Africa_Europe_IBA/Africa_Europe_IBA_centroids.shp") 
site_cent %<>% 
  mutate(SitRecID = as.character(SitRecID)) %>% 
  dplyr::select(SitRecID, geometry)

## re-merge in and out data
netdat <- outdat %>% st_drop_geometry() %>% bind_rows(indat) %>%
  arrange(bird_id, timestamp)

# create reference table of sites w/ obs --------------------------------------
site_cent <- bind_rows(site_cent, out_cent)

site_summ <- netdat %>% 
  group_by(site_poly, SitRecID) %>% summarise() %>% 
  left_join(site_cent, by = c("SitRecID")) %>% st_as_sf()


###---------------------------------------------------------------------------
### network ------------------------------------------------------------------

## create (relative) numeric code for nodes/sites -----------------------------
# netdat$loc_num <- as.numeric(as.factor(netdat$site_poly)) # name
netdat$loc_num <- as.numeric(as.factor(netdat$SitRecID))  # (absolute) numeric

## Edge list ------------------------------------------------------------------

## produce all combinations of sites visited by each individual
netdat_list <- split(netdat, netdat$bird_id)

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

nodelist <- netdat %>% group_by(loc_num, SitRecID) %>% 
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
# mapview::mapview(nodesf, zcol="degree_rank")
# mapview::mapview(nodesf, zcol="btwn_rank")
mapview::mapview(edgesf) + mapview::mapview(nodesf)

# mapview::mapview(nodesf, zcol="btwn_rank") 
# mapview::mapview(nodesf, zcol="btwn_rank") +
# (filter(edgesf, from %in% 241 | to %in% 241) %>% 
#   mapview::mapview())

## SAVE ##

# saveRDS(netsf, paste0("data/analysis/networks/", datatype,"_", season, "_iba10km_poly.rds"))
