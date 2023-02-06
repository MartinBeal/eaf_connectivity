### Map individual relocation data vs network for checking results etc.

pacman::p_load(dplyr, igraph, stringr, tictoc, tidygraph, sfnetworks, ggplot2, 
               sf, mapview, magrittr, lubridate, netrankr, dggridR, data.table)


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
if(datatype == "trax"){
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

## centroids of IBA polygons
site_cent <- readRDS( 
  "data/geodata/ibas/Africa_Europe_IBA/Africa_Europe_IBA_centroids.shp") 

if(datatype == "trax"){
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
  oneid %<>% rename(site_poly = IntName) ## IBAs only
  
  ## split data by in and outside of IBA
  outdat  <- filter(oneid, site_poly == "none")
  indat  <- filter(oneid, site_poly != "none")
  
  ## Hexgrid --------------------------------------------------------------------
  
  ## Construct global grid with cells approximately 10 km across
  # dggs <- dgconstruct(spacing=8, resround='down')
  dggs <- dgconstruct(spacing=10, resround='down')
  
  ## extract cell # for each relocation
  outdat$cell <- dgGEO_to_SEQNUM(
    dggs, outdat$longitude, outdat$latitude)$seqnum
  
  outdat$SitRecID <- paste0("cell_", as.character(outdat$cell)) # use cell # as absolute ID
  
  ## re-merge in and out data
  oneid <- indat %>% mutate(cell = NA) %>% bind_rows(outdat) %>%
    arrange(bird_id, timestamp)
  
  
  ## Calculate number of consecutive days spent at a site -----------------------
  oneid <- as.data.table(oneid)
  
  if(datatype == "trax"){
    oneid[, c("samesite","n_day"):=NULL] # remove cols
    
    ## id consecutive obs at a site/cell, per bird (data.table way)
    oneid[, samesite := data.table::rleid(SitRecID), by = bird_id]
    oneid$samesite <- ifelse(is.na(oneid$SitRecID), NA, oneid$samesite) # NA where site NA
    
    ## number of consecutive days at a site
    oneid <- oneid %>%
      group_by(bird_id, samesite) %>%
      summarise(n_day = n_distinct(yday(timestamp))) %>%
      left_join(oneid)
    
    ## for tracking data, remove data from sites visited for < 48h
    oneid <- filter(oneid, n_day >= 2)
  }
  
  ## check out some data
  # oneid %>% filter(bird_id == "B.tenskar") %>% 
  # oneid[5000:10000,] %>%
  #   sf::st_as_sf(
  #     coords = c("longitude", "latitude"),
  #     crs = 4326, agr = "constant") %>% mapview(zcol="SitRecID")
  
  ## split data again by in and out site to add cell center coords --------------
  outdat <- oneid %>% 
    filter(site_poly == "none")
  indat <- oneid %>% 
    filter(site_poly != "none")
  
  ## SKip if ind has no points outside IBAs -----------------------------------
  if(nrow(outdat) > 0){
    ## one row per cell
    binned_re <- outdat %>% group_by(SitRecID, cell) %>% 
      summarise() # %>% summarise(
    # count=n(), 
    # n_birds = n_distinct(bird_id))
    
    ## Get the grid cell boundaries for cells which had obs
    grid_re           <- dgcellstogrid(dggs, binned_re$cell)
    colnames(grid_re)[1] <- "cell"
    
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
        site_type = "outsite",
        cellgrp = paste0("cellgrp_", cellgrp),
        SitRecID = cellgrp ## make the cell group the main site ID
      ) %>% group_by(SitRecID, cellgrp, site_type) %>% 
      st_as_sf(coords = c("longitude", "latitude"), 
               crs = 4326, agr = "constant") %>% 
      summarise() %>% 
      dplyr::select(SitRecID, site_type, geometry)
    
    site_cent %<>% 
      mutate(
        site_type = "IBA",
        SitRecID = as.character(SitRecID)) %>% 
      dplyr::select(SitRecID, site_type, geometry)
    
    ## re-merge in and out data
    oneid <- outdat %>% 
      mutate(longitude = st_coordinates(.)[,1],
             latitude = st_coordinates(.)[,2]) %>% 
      st_drop_geometry() %>% bind_rows(indat) %>%
      arrange(bird_id, timestamp)
    
    # create reference table of sites w/ obs --------------------------------------
    site_cent <- bind_rows(site_cent, out_cent)
  } else {
    site_cent %<>% 
      mutate(
        site_type = "IBA",
        SitRecID = as.character(SitRecID)) %>% 
      dplyr::select(SitRecID, site_type, geometry)
    
    ## re-merge in and out data
    oneid <- indat
  }

  site_summ <- oneid %>% 
    group_by(site_poly, SitRecID) %>% summarise() %>% 
    left_join(site_cent, by = c("SitRecID")) %>% st_as_sf()
  
  site_summ %>% mapview(zcol="site_type")
  
  if(nrow(site_summ) < 2){
    print("only one site visited")
    next}
  
  ## Edge list ------------------------------------------------------------------
  
  ###---------------------------------------------------------------------------
  ### network ------------------------------------------------------------------
  
  ## create (relative) numeric code for nodes/sites -----------------------------
  # oneid$loc_num <- as.numeric(as.factor(oneid$site_poly)) # name
  oneid$loc_num <- as.numeric(as.factor(oneid$SitRecID))  # (absolute) numeric
  
  
  ## Edge list ------------------------------------------------------------------
  
  ## produce all combinations of sites visited by each individual
  oneid_list <- split(oneid, oneid$bird_id)
  
  oneid_list <- lapply(
    seq_along(oneid_list), 
    function(x){
      # print(x)
      one <- oneid_list[[x]]
      xx <- as.data.frame(
        RcppAlgos::comboGrid(
          one$loc_num,
          one$loc_num, repetition = TRUE)
      )
      xx$bird_id <- one$bird_id[1]
      return(xx)
    })
  
  full <- data.table::rbindlist(oneid_list)
  
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
    n_obstype <- oneid %>% 
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
    n_obstype <- oneid %>% 
      group_by(loc_num, SitRecID, device) %>% 
      summarise(
        n_obs = n()
      ) %>% tidyr::pivot_wider(
        names_from = "device",
        names_prefix = "n_",
        values_from = "n_obs"
      )
  }
  
  n_id_total <- n_distinct(oneid$bird_id)
  
  nodelist <- oneid %>% group_by(loc_num, SitRecID) %>% 
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
  
  
  oneid_sf <- st_as_sf(oneid, coords = c("longitude", "latitude"), 
                       crs = 4326, agr = "constant")
  
  nodesf <- netsf %>% activate("nodes") %>% sf::st_as_sf()
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
    geom_sf(data = nodesf, aes(), color="red") + 
    geom_sf(data = oneid_sf, aes(), color="blue", size = .35, alpha=.5) + 
    coord_sf(xlim = c(bbox[1], bbox[3]), 
             ylim = c(bbox[2], bbox[4]), expand = T) + 
    theme_void() +
    theme(panel.background = element_rect(fill="white"))# +
  # ggtitle("Network")
  
  id <- oneid$bird_id[1]
  
  # library(patchwork)
  # combmap <- trxmap + netmap
  # 
  # ggsave(paste0("figures/ind_data/", datatype, "/", season, "/", id, ".png"), 
  #        height=3.5, width=3.5)
  # 
  ggsave(paste0("figures/ind_data/", datatype, "/", season, "/", id, ".png"),
         plot=netmap, height=3.5, width=3.5)
}

