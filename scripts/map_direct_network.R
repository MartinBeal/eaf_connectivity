## Map DIRECTED, SEASONAL networks ---------------------------------

pacman::p_load(sfnetworks, sf, dplyr, magrittr, ggplot2, stringr, units)

## which network to create
# season <- "all"
# season <- "spring"
season <- "fall"

## networks including in and out sites
anet <- readRDS(
  paste0("data/analysis/networks/directed/32d_", season, "_iba_hex_10km.rds")
  )

## interactive map it
nodesf <- anet %>% activate("nodes") %>% sf::st_as_sf()
edgesf <- anet %>% activate("edges") %>% sf::st_as_sf()
# mapview::mapview(nodesf, zcol="n_id")
# mapview::mapview(edgesf)

## Dataframe of PA coverage of sites/nodes
site_cover <- readRDS("data/analysis/protectedness/perc_pa_cover_allsites.rds")

## add in protectedness info
nodesf <- nodesf %>% left_join(site_cover)

## summarise seasonal coverage across flyway ----------------------------------
season_summ <- nodesf %>% st_drop_geometry() %>% 
  summarise(
    n_sites     = n(),
    n_anycover  = sum(ifelse(as.numeric(na.omit(pct_cover)) > 0, T, F)),
    perc_wcover = (n_anycover / n_sites)*100,
    mn_cover = mean(na.omit(pct_cover)),
    sd_cover = sd(na.omit(pct_cover)),
  )

## SAVE

saveRDS(
  season_summ, paste0("data/analysis/protectedness/perc_pa_cover_", season, ".rds")
  )


## Maps: project for prettier map --------------------------------------

## site polygons
# rams <- raster::shapefile("data/geodata/ramsar/EAF_ramsar.shp")
# iba <- raster::shapefile("data/geodata/ibas/EAF_btgo_IBA/eaf_btgo_iba.shp")

# get world map
wmap <- rworldmap::getMap(resolution="high")
wmap_prj <- st_as_sf(wmap) %>% st_transform(crs = "EPSG:3035")

## set common theme to re-use
maptheme <- theme(
  plot.background = element_rect(fill = "white"),
  panel.background = element_rect(fill = "white"),
  axis.text = element_blank(),
  axis.ticks = element_blank(),
  panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
  panel.border = element_rect(colour = "black", fill=NA, size=.5),
  legend.key=element_blank()
)

## Lambert EA: EPSG:3035 ## Lambert CConic: EPSG:3034
node_prj <- nodesf %>% st_transform("EPSG:3035")
edge_prj <- edgesf %>% st_transform("EPSG:3035")

## EAF bbox
bbox_prj <- st_bbox(
  c(xmin = -16, xmax = 21,
    ymin = 8, ymax = 65.2), crs = 4326) %>% 
  st_as_sfc() %>% st_transform("EPSG:3035") %>% st_bbox()

## deal w/ NAs
# edge_prj$n_bird_trax <- ifelse(
#   is.na(edge_prj$n_bird_trax), -1, edge_prj$n_bird_trax)
  
## color edges IDd including tracking data and those from only ring data
edge_prj$dtype_grps <- ifelse(
  str_detect(edge_prj$datatypes, pattern="trax"), "Tracking", "Ring data")
  
## ID disconnected nodes
disc_nod_ids <- unique(node_prj$SitRecID)[
  which(!unique(node_prj$SitRecID) %in% as.character(
    unique(c(edge_prj$from_sid, edge_prj$to_sid))
    ))]
node_prj$connected <- ifelse(
  node_prj$SitRecID %in% disc_nod_ids, "No", "Yes")
  

  ### Map 1a: edges colored by data origin, nodes by disconnectedness ---------
  map1a <- ggplot() +
    geom_sf(data = wmap_prj, fill = "grey70", color = NA) +
    # geom_sf(data = edge_prj,                 ## edges
    #         aes(size = n_bird_trax), col = "black") + #, alpha = 0.65
    geom_sf(data = edge_prj,                 ## edges
            aes(color = dtype_grps)) + #, alpha = 0.65
    geom_sf(data = wmap_prj, fill = NA, color = "white", size=0.2) +
    geom_sf(data = arrange(node_prj, n_id), 
            aes(fill = connected),
            pch=21, stroke=.25,
            col="black") +
    coord_sf(xlim = 
               c(bbox_prj[1], bbox_prj[3]), ylim = c(bbox_prj[2], bbox_prj[4]), 
             expand = T) +
    scale_fill_brewer(type="qual") +
    labs(color = "Data type",
         fill = "Connected") +
  # scale_size(range = c(0.01, 2)) +
    maptheme
  map1a
  
  ## SAVE  ##
  ggsave(
    paste0("figures/direct_networks/dtypeedge_discnode_", season, ".png"), 
    plot=map1a, width=5, height = 6)

  
  ### Map 1b: edges colored by data origin, rmv few ID disconnect. nodes ------
  node_prj2 <- filter(node_prj, (connected == "No" & n_id > 1) |
                        (connected == "Yes"))
  
  map1b <- ggplot() +
    geom_sf(data = wmap_prj, fill = "grey70", color = NA) +
    # geom_sf(data = edge_prj,                 ## edges
    #         aes(size = n_bird_trax), col = "black") + #, alpha = 0.65
    geom_sf(data = edge_prj,                 ## edges
            aes(color = dtype_grps)) + #, alpha = 0.65
    geom_sf(data = wmap_prj, fill = NA, color = "white", size=0.2) +
    geom_sf(data = arrange(node_prj2, n_id), 
            aes(fill = connected),
            pch=21, stroke=.25,
            col="black") +
    # geom_sf(data = arrange(node_prj2, n_id), size=1) +
    coord_sf(xlim = 
               c(bbox_prj[1], bbox_prj[3]), ylim = c(bbox_prj[2], bbox_prj[4]), 
             expand = T) +
    scale_fill_brewer(type="qual") +
    # scale_size(range = c(0.01, 2)) +
    labs(color = "Data type",
         fill = "Connected") +
    maptheme
  map1b
  
  ## SAVE  ##
  ggsave(
    paste0("figures/direct_networks/dtypeedge_fewdiscnode_", season, ".png"), 
    plot=map1b, width=5, height = 6)
  
  
  ### Map 1c: edges colored by data origin, rmv few ID disconnect. nodes ------
  map1c <- ggplot() +
    geom_sf(data = wmap_prj, fill = "grey70", color = NA) +
    geom_sf(data = edge_prj) + #, alpha = 0.65
    geom_sf(data = wmap_prj, fill = NA, color = "white", size=0.2) +
    # geom_sf(data = node_prj2, pch=21, fill="white", col="black") +
    coord_sf(xlim = 
               c(bbox_prj[1], bbox_prj[3]), ylim = c(bbox_prj[2], bbox_prj[4]), 
             expand = T) +
    scale_fill_brewer(type="qual") +
    maptheme
  map1c
  
  ## SAVE  ##
  ggsave(
    paste0("figures/direct_networks/onlyedges_", season, ".png"), 
    plot=map1c, width=5, height = 6)
  
  
  ## Map 2:3 network metrics --------------------------------------------------
  ## top 25 sites - IN DEGREE 
  node_top_degin <- arrange(node_prj, desc(degree_in)) %>% slice(1:25)
  
  map2a <- ggplot() +
    geom_sf(data = wmap_prj, fill = "grey70", color = NA) +
    # geom_sf(data = iba_sf_prj,               ## IBAs
    #         aes(), fill = "dark red", color = NA) +
    # geom_sf(data = edge_prj,                 ## edges
    #         aes(size = degree), col = "black") + #, alpha = 0.65
    geom_sf(data = wmap_prj, fill = NA, color = "white", size=0.2) +
    # geom_sf(data = arrange(node_prj, between), ## node betweeness
    #         aes(color = between, size = between)) +
    # geom_sf(data = arrange(node_20_deg, degree),  ## node degree
    #         aes(color = degree, size = degree)) +
    geom_sf(data = arrange(node_top_degin, degree_in_rank),  ## node degree rank
            aes(color = degree_in_rank, size = degree_in_rank)) +
    coord_sf(
      xlim = c(bbox_prj[1], bbox_prj[3]), ylim = c(bbox_prj[2], bbox_prj[4]), 
      expand = T) +
    scale_size(range = c(0.01, 2),   # reverse size scale (for ranks)
               trans = 'reverse') +
    maptheme
  # reverse size scale (for ranks)
  # scale_size(range = c(0.01, 2)) 
  map2a
  
  ## top 25 sites - OUT DEGREE 
  node_top_degout <- arrange(node_prj, desc(degree_out)) %>% slice(1:25)
  
  map2b <- ggplot() +
    geom_sf(data = wmap_prj, fill = "grey70", color = NA) +
    geom_sf(data = wmap_prj, fill = NA, color = "white", size=0.2) +
    geom_sf(data = arrange(node_top_degout, degree_out_rank),  ## node degree rank
            aes(color = degree_out_rank, size = degree_out_rank)) +
    coord_sf(
      xlim = c(bbox_prj[1], bbox_prj[3]), ylim = c(bbox_prj[2], bbox_prj[4]), 
      expand = T) +
    scale_size(range = c(0.01, 2),   # reverse size scale (for ranks)
               trans = 'reverse') +
    maptheme
  # reverse size scale (for ranks)
  # scale_size(range = c(0.01, 2)) 
  map2b
  
  ## top 25 sites - OUT DEGREE 
  node_top_degtot <- arrange(node_prj, desc(degree_tot)) %>% slice(1:25)
  
  map2c <- ggplot() +
    geom_sf(data = wmap_prj, fill = "grey70", color = NA) +
    geom_sf(data = wmap_prj, fill = NA, color = "white", size=0.2) +
    geom_sf(data = arrange(node_top_degtot, degree_tot_rank),  ## node degree rank
            aes(color = degree_tot_rank, size = degree_tot_rank)) +
    coord_sf(
      xlim = c(bbox_prj[1], bbox_prj[3]), ylim = c(bbox_prj[2], bbox_prj[4]), 
      expand = T) +
    scale_size(range = c(0.01, 2),   # reverse size scale (for ranks)
               trans = 'reverse') +
    maptheme
  # reverse size scale (for ranks)
  # scale_size(range = c(0.01, 2)) 
  map2c
  
  ## top 25 sites - BETWEENNESS 
  node_top_btwn <- arrange(node_prj, desc(between)) %>% slice(1:25)
  
  map3 <- ggplot() +
    geom_sf(data = wmap_prj, fill = "grey70", color = NA) +
    # geom_sf(data = iba_sf_prj,               ## IBAs
    #         aes(), fill = "dark red", color = NA) +
    # geom_sf(data = edge_prj,                 ## edges
    #         aes(size = degree), col = "black") + #, alpha = 0.65
    geom_sf(data = wmap_prj, fill = NA, color = "white", size=0.2) +
    # geom_sf(data = arrange(node_prj, between), ## node betweeness
    #         aes(color = between, size = between)) +
    # geom_sf(data = arrange(node_prj, degree),  ## node degree
    #         aes(color = degree, size = degree)) +
    geom_sf(data = arrange(node_top_btwn, between),  ## node degree rank
            aes(color = btwn_rank, size = btwn_rank)) +
    coord_sf(
      xlim = c(bbox_prj[1], bbox_prj[3]), ylim = c(bbox_prj[2], bbox_prj[4]), 
      expand = T) +
    scale_size(range = c(0.01, 2), # reverse size scale (for ranks) 
               trans = 'reverse') +
    maptheme
  # scale_size(range = c(0.01, 2)) +
  map3
  
  ## Save -------------------------------------
  ggsave(
    paste0("figures/direct_networks/", season, "_tot_degreernkX.png"), 
    plot=map2a, width=5, height = 6)
  ggsave(
    paste0("figures/direct_networks/", season, "_in_degreernkX.png"), 
    plot=map2b, width=5, height = 6)
  ggsave(
    paste0("figures/direct_networks/", season, "_out_degreernkX.png"), 
    plot=map2c, width=5, height = 6)
  ggsave(
    paste0("figures/direct_networks/", season, "_btwnrnkX.png"), 
    plot=map3, width=5, height = 6)
  
  
  ### Protectedness 
  map4 <- ggplot() +
    geom_sf(data = wmap_prj, fill = "grey70", color = NA) +
    geom_sf(data = wmap_prj, fill = NA, color = "white", size=0.2) +
    geom_sf(data = arrange(node_prj, degree),  ## node degree rank
            aes(fill = pct_cover), pch=21) +
    scale_fill_distiller(
      palette="Reds",
      guide = guide_legend(
        override.aes = list(size=4) ## make color legend points not bar
      )) +
    coord_sf(
      xlim = c(bbox_prj[1], bbox_prj[3]), ylim = c(bbox_prj[2], bbox_prj[4]), 
      expand = T) +
    scale_size(range = c(0.01, 2),   # reverse size scale (for ranks)
               trans = 'reverse') +
    maptheme 
  # reverse size scale (for ranks)
  # scale_size(range = c(0.01, 2))
  map4
  
  ## SAVE
  ggsave(
    paste0("figures/direct_networks/", season, "_pa_cover.png"), 
    plot=map4, width=5, height = 6)
  