### Create a spatial, undirected network from ring recapture and resighting data
## to sites defined using the IBA layer and hexcells for unidentified sites

# VERSION that uses same sites (in and out) for all data/seasons

pacman::p_load(dplyr, igraph, stringr, tictoc, tidygraph, sfnetworks, ggplot2, 
               sf, mapview, magrittr, lubridate, netrankr)


## choose data subset to run --------------------------------------------------

## which network to create
# season <- "all"
# season <- "spring"
season <- "fall"

## Run through each data type -------------------------------------------------
# dtype <- "metal"
# dtype <- "color"
dtype <- "trax"

## Load site cummary both IBAs and outsite centroids) ------------------------
site_summ <- readRDS(
  paste0("data/analysis/site_nodes/alldatatypes_allsites_cent_all.rds")
)

## Load combined data, with both IBA overlay info and hexgrid outsites -------
### SAVE 
alldat <- readRDS(
  paste0("data/analysis/combined/alldatatypes_ibas_outsites_all.rds"))

## filter to one datatype
alldat <- subset(alldat, str_detect(alldat$datatype, dtype))

# fxn for splitting string into columns
source("C:/Users/Martim Bill/Documents/R/source_scripts/str2col.R")
# fxn for calc. displacement distance from first tracking locations
source("C:/Users/Martim Bill/Documents/R/source_scripts/displ_dist.R")

# alldat %<>% rename(combo = metal)
# alldat %<>% filter(obssource == "euring") # which data sources
# alldat %<>% filter(obssource != "euring")

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
if(dtype == "trax"){
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
  
  wmig <- filter(displ_id, mx_displ >= 100)
  
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


###---------------------------------------------------------------------------
### network ------------------------------------------------------------------

## remove birds w/ only one sighting (to avoid unconnected nodes) -------
xz <- netdat %>% group_by(bird_id) %>% 
  summarise(nobs = n()) %>% filter(nobs == 1)
netdat <- subset(netdat, !bird_id %in% xz$bird_id)

## remove birds only seen (multiple times) at same site ----------------------- 
nsites <- netdat %>% group_by(bird_id) %>% 
  summarise(nsites = n_distinct(site_poly))

## % of birds w/ relocs at one site only
sum(nsites$nsites == 1) / n_distinct(netdat$bird_id) * 100

xz2 <- filter(nsites, nsites == 1)
netdat <- subset(netdat, !bird_id %in% xz2$bird_id)

## (local) numeric code of nodes/sites ------------------------------------------------
# netdat$loc_num <- as.numeric(as.factor(netdat$site_poly)) # name
netdat$loc_num <- as.numeric(as.factor(netdat$SitRecID))  # (absolute) numeric

## Edge list ------------------------------------------------------------------

oneid_list <- split(netdat, netdat$bird_id)

oneid_list <- lapply(
  seq_along(oneid_list), 
  function(x){
    one <- oneid_list[[x]]
    xx <- data.frame(
      bird_id = one$bird_id[i],
      from = one$loc_num[1:nrow(one)-1],
      to   = one$loc_num[2:nrow(one)]
    )
    return(xx)
  })

full <- data.table::rbindlist(oneid_list)

## remove self connections
noself <- full[-which(full$from == full$to), ]

## combine sites into single variable for summarizing
noself$sitecomb <- paste(noself$from, noself$to)

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
  dplyr::select(from, to, n_id, prop_id) %>% 
  as_tibble()


## Vertex list ---------------------------------------------------------------
# if(dtype %in% c("color", "metal")){
#   n_obstype <- netdat %>% 
#     group_by(loc_num, SitRecID, obstype) %>% 
#     summarise(
#       n_obs = n()
#     ) %>% tidyr::pivot_wider(
#       # id_cols = "obstype",
#       names_from = "obstype",
#       names_prefix = "n_",
#       values_from = "n_obs"
#     )
# } else if (dtype == "trax"){
#   n_obstype <- netdat %>% 
#     group_by(loc_num, SitRecID, device) %>% 
#     summarise(
#       n_obs = n()
#     ) %>% tidyr::pivot_wider(
#       names_from = "device",
#       names_prefix = "n_",
#       values_from = "n_obs"
#     )
# }

n_id_total <- n_distinct(netdat$bird_id)

nodelist <- netdat %>% group_by(loc_num, SitRecID) %>% 
  summarise(
    n_id    = n_distinct(bird_id),
    prop_id = n_id / n_id_total,
    n_obs   = n()
  ) #%>% right_join(n_obstype)

nodelist <- nodelist %>% 
  left_join(site_summ) %>% 
  sf::st_as_sf()

## Convert to sfnetwork -------------------------------------------------------

## undirected - reciprocal edges
netsf <- sfnetwork(nodelist, edgelist, node_key = "dist", directed = F, 
                   edges_as_lines = TRUE)
## directed - order matters
# netsf <- sfnetwork(nodelist, edgelist, node_key = "dist", directed = T, 
#                    edges_as_lines = TRUE)


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

saveRDS(netsf, paste0("data/analysis/networks/", dtype,"_", season, "_iba_hex_10km.rds"))

