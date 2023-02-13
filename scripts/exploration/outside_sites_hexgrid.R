#### Identify used sites/nodes and their importance from data falling outside
### polygon  dataset of known sites
## Sites defined using hexgrid binning (i.e. coordinates binned into cells
# of a certain size)

pacman::p_load(dplyr, igraph, stringr, tictoc, tidygraph, sfnetworks, ggplot2,
               sf, mapview, magrittr, lubridate, dggridR, data.table)

### Data type
# datatype <- "metal"
# datatype <- "color"
datatype <- "trax"

### Season
season <- "all"
# season <- "spring"
# season <- "fall"

if(datatype == "color"){
  ## ring relocations overlaid on polygon layer
  alldat <- readRDS(
    paste0("data/analysis/ringing/color_outside_", season, "_ibas10km.rds"))
} else if (datatype == "metal"){
  ## metal ring locations overlaid on polygon layer
  alldat <- readRDS(
    paste0("data/analysis/ringing/metal_outside_", season, "_ibas10km.rds"))
} else if (datatype == "trax"){
  ## tracking locations overlaid on polygon layer
  alldat <- readRDS(
    paste0("data/analysis/tracking/outside_trx_", season, "_ibas10km.rds"))
}

## fxn for splitting string into columns
source("C:/Users/Martim Bill/Documents/R/source_scripts/str2col.R")


## Hexgrid --------------------------------------------------------------------

netdat <- alldat ## all individuals

## Construct global grid with cells approximately 1000 miles across
# dggs <- dgconstruct(spacing=8, resround='down')
dggs <- dgconstruct(spacing=10, resround='down')

## create hexgrid from cr obs
netdat$cell <- dgGEO_to_SEQNUM(
  dggs, netdat$longitude, netdat$latitude)$seqnum

#Get the number of records in each equally-sized cell
binned_re <- netdat %>% group_by(cell) %>% summarise(count=n())

#Get the grid cell boundaries for cells which had obs
grid_re           <- dgcellstogrid(dggs, binned_re$cell)
colnames(grid_re)[1] <- "cell"
ccenters_re       <- dgSEQNUM_to_GEO(dggs, grid_re$cell)

# cbind.data.frame(lat=ccenters_re$lat_deg, lon=ccenters_re$lon_deg) %>%
#   sf::st_as_sf(
#     coords = c("lon", "lat"),
#     crs = 4326, agr = "constant") %>% mapview()

#Update the grid cells' properties to include the number of obs in each cell
grid_re <- merge(grid_re, binned_re, by.x="cell", by.y="cell")

# mapview(grid_re, zcol = "count")

## add cell center coordinates
grid_re <- grid_re %>% st_drop_geometry() %>%
  bind_cols(latitude=ccenters_re$lat_deg, longitude=ccenters_re$lon_deg)


## Calculate number of consecutive days spent at a site -----------------------
netdat <- as.data.table(netdat)
if(datatype == "trax"){
  netdat[, c("samesite","n_day"):=NULL]
}

## id consecutive obs at a site, per bird (data.table way)
netdat[, samesite := data.table::rleid(cell),
       by = bird_id]
netdat$samesite <- ifelse(is.na(netdat$cell), NA, netdat$samesite) # NA where site NA

## number of consecutive days at a site
netdat <- netdat %>%
  group_by(bird_id, samesite) %>%
  summarise(n_day = n_distinct(yday(timestamp))) %>%
  left_join(netdat)

## for tracking data, data outside cells visited for < 48h
if (datatype == "trax"){
  netdat <- filter(netdat, n_day >= 2)
}


###---------------------------------------------------------------------------
### network ------------------------------------------------------------------

## numeric coding of locations
netdat$loc_num <- as.numeric(as.factor(netdat$cell))

netdat <- netdat %>%
  dplyr::select(-latitude, -longitude) %>%
  left_join(grid_re, by="cell")

# create reference table of sites w/ obs
site_summ <- netdat %>% group_by(loc_num) %>%
  summarise(latitude = first(latitude), longitude = first(longitude))


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
if(datatype %in% c("color", "metal")){
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
} else if (datatype == "trax"){
  n_obstype <- netdat %>%
    group_by(loc_num, device) %>%
    summarise(
      n_obs = n()
    ) %>% tidyr::pivot_wider(
      names_from = "device",
      names_prefix = "n_",
      values_from = "n_obs"
    )
}

nodelist <- netdat %>% group_by(loc_num) %>%
  summarise(
    n_id = n_distinct(bird_id),
    n_obs  = n()
  ) %>% right_join(n_obstype)

nodelist <- nodelist %>%
  left_join(site_summ, by="loc_num") %>%
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
netsf <- netsf %>%
  activate("nodes") %>%
  filter(n_id > 1)

## interactive map it
nodesf <- netsf %>% activate("nodes") %>% sf::st_as_sf()
# netsf %>% activate("edges") %>% sf::st_as_sf() %>% mapview::mapview() +
#   mapview::mapview(nodesf, zcol="n_id")
## just nodes
mapview::mapview(nodesf, zcol="n_id")


## SAVE ## 

if(datatype == "color"){
  saveRDS(netsf, 
          paste0("data/analysis/networks/color_outside_", season, "_iba10km_poly.rds"))
} else if (datatype == "metal"){
  saveRDS(netsf, 
          paste0("data/analysis/networks/metal_outside_", season, "_iba10km_poly.rds"))
} else if (datatype == "trax"){
  saveRDS(netsf, 
          paste0("data/analysis/networks/trax_outside_", season, "_iba10km_poly.rds"))
}

