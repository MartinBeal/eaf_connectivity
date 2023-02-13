### Map individual relocation data vs network for checking results etc.

pacman::p_load(dplyr, igraph, stringr, tictoc, tidygraph, sfnetworks, ggplot2, 
               sf, mapview, magrittr, lubridate, netrankr, dggridR, data.table)


## Run through each data type -------------------------------------------------
# dtype <- "metal"
# dtype <- "color"
dtype <- "trax"

# if(dtype == "color"){
#   ## color ringed bird captures and resightings overlaid on polygon layer
#   alldat <- readRDS("data/analysis/ringing/cring_merge_no7dayreobs_ibas.rds")
#   alldat %<>% rename(timestamp = date)
# } else if (dtype == "metal"){
#   ## metal ring captures, recaptures, recoveries
#   alldat <- readRDS("data/analysis/ringing/euring_metal_ibas.rds")
#   alldat %<>% rename(timestamp = date)
# } else if (dtype == "trax"){
#   ## tracking locations 
#   # alldat <- readRDS("data/analysis/tracking/PTT_GPS_mconn_12h_ibas.rds")
#   alldat <- readRDS("data/analysis/tracking/PTT_GPS_mconn_12h_no0_ibas.rds")
#   # alldat %<>% rename(timestamp = date)
# }

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


## choose data subset to run --------------------------------------------------

## which network to create
season <- "all"
# season <- "spring"
# season <- "fall"

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
  
  ## Summarise max/avg distance from first points -----------------------------
  
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

#### Run one individual only --------------------------------------------------

n_ids <- n_distinct(netdat$bird_id)

## centroids of all sites (in and out)
site_summ <- readRDS(
  "data/analysis/site_nodes/alldatatypes_allsites_cent_all.rds"
)

if(dtype == "trax"){
  idxs <- seq_len(n_ids)
} else {
  idxs <- 1:300
}
# for(i in seq_len(n_ids)){
for(i in idxs){
  
  print(i)
  
  oneid <- filter(netdat, bird_id == unique(netdat$bird_id)[i])
  
  ## top sites visited
  # xx <- oneid2 %>% group_by(SitRecID, IntName, country) %>% 
  #   summarise(
  #     n_birds = n_distinct(bird_id),
  #     n_obs   = n()
  #   ) %>% arrange(desc(n_birds)) %>% ungroup() %>% 
  #   filter(country != "Iceland") %>% 
  #   slice(1:11)
  
  ## 
  # oneid %<>% rename(site_poly = IntName) ## IBAs only
  
  ## split data by in and outside of IBA
  outdat  <- filter(oneid, site_poly == "none")
  indat  <- filter(oneid, site_poly != "none")
  
  ## check out some data
  # oneid %>% filter(bird_id == "B.tenskar") %>% 
  # oneid[5000:10000,] %>%
  #   sf::st_as_sf(
  #     coords = c("longitude", "latitude"),
  #     crs = 4326, agr = "constant") %>% mapview(zcol="SitRecID")
  
  if(n_distinct(oneid$SitRecID) < 2){
    print("only one site visited")
    next}
  
  ## Edge list ------------------------------------------------------------------
  
  ###---------------------------------------------------------------------------
  ### network ------------------------------------------------------------------
  
  ## create (relative) numeric code for nodes/sites -----------------------------
  # oneid$loc_num <- as.numeric(as.factor(oneid$site_poly)) # name
  oneid$loc_num <- as.numeric(as.factor(oneid$SitRecID))  # (absolute) numeric
  
  
  ## Edge list ------------------------------------------------------------------
  oneid_list <- split(oneid, oneid$bird_id)
  
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
  #   n_obstype <- oneid %>% 
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
  #   n_obstype <- oneid %>% 
  #     group_by(loc_num, SitRecID, device) %>% 
  #     summarise(
  #       n_obs = n()
  #     ) %>% tidyr::pivot_wider(
  #       names_from = "device",
  #       names_prefix = "n_",
  #       values_from = "n_obs"
  #     )
  # }
  
  n_id_total <- n_distinct(oneid$bird_id)
  
  nodelist <- oneid %>% group_by(loc_num, SitRecID) %>% 
    summarise(
      n_id    = n_distinct(bird_id),
      prop_id = n_id / n_id_total,
      n_obs   = n()
    ) #%>% right_join(n_obstype)
  
  nodelist <- nodelist %>% 
    left_join(site_summ) %>% 
    sf::st_as_sf()
  
  
  ## Convert to sfnetwork -------------------------------------------------------
  
  # netsf <- sfnetwork(nodelist, edgelist, node_key = "dist", directed = F)
  netsf <- sfnetwork(nodelist, edgelist, node_key = "dist", directed = F, 
                     edges_as_lines = TRUE)
  
  
  oneid_sf <- st_as_sf(oneid, coords = c("longitude", "latitude"), 
                       crs = 4326, agr = "constant")
  
  nodesf <- netsf %>% activate("nodes") %>% sf::st_as_sf() %>% mutate(
    site_type = ifelse(site_poly == "none", "Outsite", "IBA")
  )
  edgesf <- netsf %>% activate("edges") %>% sf::st_as_sf() %>% st_make_valid()
  
  ## 
  # get world map
  wmap <- rworldmap::getMap(resolution="low")
  wmap_prj <- st_as_sf(wmap) 
  
  bbox <- st_bbox(
    c(xmin = -23, xmax = 30,
      ymin = 7, ymax = 66), crs = 4326)
  
  # trxmap <- ggplot() + 
  #   geom_sf(data = wmap_prj, fill = NA, color = "grey20", size=0.2) +
  #   geom_sf(data = oneid_sf, aes(), color="red", alpha=.5) + 
  #   coord_sf(xlim = c(bbox[1], bbox[3]), 
  #            ylim = c(bbox[2], bbox[4]), expand = T) + 
  #   theme_void() +
  #   theme(panel.background = element_rect(fill="white")) +
  #   
  #   ggtitle("Relocations")
  
  netmap <- ggplot() + 
    geom_sf(data = wmap_prj, fill = NA, color = "grey20", size=0.2) +
    geom_sf(data = edgesf, aes()) +
    # geom_sf(data = nodesf, aes(), color="red") + 
    geom_sf(data = nodesf, aes(color=site_type)) + 
    geom_sf(data = oneid_sf, aes(), color="purple", size = .35, alpha=.5) + 
    coord_sf(xlim = c(bbox[1], bbox[3]), 
             ylim = c(bbox[2], bbox[4]), expand = T) + 
    theme_void() +
    theme(
      panel.background = element_rect(fill="white"),
      legend.background = element_rect(fill="white", color = NA),
      legend.key = element_rect(fill="white", color = NA))# +
  # ggtitle("Network")
  
  id <- oneid$bird_id[1]
  
  # library(patchwork)
  # combmap <- trxmap + netmap
  # 
  # ggsave(paste0("figures/ind_data/", dtype, "/", season, "/", id, ".png"), 
  #        height=3.5, width=3.5)
  # 
  ggsave(paste0("figures/ind_data/", dtype, "/", season, "/", id, ".png"),
         plot=netmap, height=3.5, width=3.5)
}

