## Compare rank differences between networks of differing data types

pacman::p_load(sfnetworks, sf, dplyr, magrittr, mapview, ggplot2, patchwork, stringr)


## loop through each season category
seasons <- c("all", "spring", "fall")

for(i in seq_along(seasons)){
  print(i)
  
  ### which network to create (manual)
  # season <- "all"
  # season <- "spring"
  # season <- "fall"
  
  season <- seasons[i]
  
  m_net <- readRDS(
    paste0("data/analysis/networks/metal_", season, "_iba_hex_10km.rds"))
  c_net <- readRDS(
    paste0("data/analysis/networks/color_", season, "_iba_hex_10km.rds"))
  t_net <- readRDS(
    paste0("data/analysis/networks/trax_", season, "_iba_hex_10km.rds"))
  
  
  ## interactive map it
  anet <- m_net
  nodesf <- anet %>% activate("nodes") %>% sf::st_as_sf()
  edgesf <- anet %>% activate("edges") %>% sf::st_as_sf()
  mapview::mapview(nodesf, zcol="n_id")
  mapview::mapview(edgesf)
  
  ## split up nodes and edges ---------------------------------------------------
  m_nodes <- st_as_sf(m_net, "nodes")
  c_nodes <- st_as_sf(c_net, "nodes")
  t_nodes <- st_as_sf(t_net, "nodes")
  
  m_edges <- st_as_sf(m_net, "edges")
  c_edges <- st_as_sf(c_net, "edges")
  t_edges <- st_as_sf(t_net, "edges")
  
  global_compare <- data.frame(
    data_type = c("together", "metal", "color", "tracking"),
    season = rep(season),
    n_nodes = c(
      n_distinct(c(m_nodes$loc_num, c_nodes$loc_num, t_nodes$loc_num)),
      n_distinct(m_nodes$loc_num), 
      n_distinct(c_nodes$loc_num), 
      n_distinct(t_nodes$loc_num)),
    n_edges = 
      c(n_distinct(c(m_edges$link_id, c_edges$link_id, t_edges$link_id)),
        n_distinct(m_edges$link_id), 
        n_distinct(c_edges$link_id), 
        n_distinct(t_edges$link_id))
  )
  global_compare
  
  ## SAVE
  write.csv(
    global_compare, 
    paste0("data/analysis/summaries/global_net_compare_", season, ".csv"), 
    row.names = F)
  
  
  ## Correspondence matrix of sites used by data types --------------------------
  comp_df <- expand.grid(
    c("metal", "color", "tracking"), c("metal", "color", "tracking"))
  
  comp_df %<>% 
    dplyr::select(Var2, Var1) %>% 
    rename(focal = Var2, coveredby=Var1) %>% 
    mutate(season = season)
  
  ## how well does data B cover data A? 
  # eg. how many sites identified from tracking are id'd from color resights?
  comp_df$n_sites_cover <- c(
    n_distinct(m_nodes$loc_num),
    sum(c_nodes$loc_num %in% m_nodes$loc_num),
    sum(t_nodes$loc_num %in% m_nodes$loc_num),
    sum(m_nodes$loc_num %in% c_nodes$loc_num),
    n_distinct(c_nodes$loc_num),
    sum(t_nodes$loc_num %in% c_nodes$loc_num),
    sum(m_nodes$loc_num %in% t_nodes$loc_num),
    sum(c_nodes$loc_num %in% t_nodes$loc_num),
    n_distinct(t_nodes$loc_num)
  )
  
  comp_df %<>% mutate(
    n_sites = ifelse(focal == "metal", n_distinct(m_nodes$loc_num),
                     ifelse(focal == "color", n_distinct(c_nodes$loc_num),
                            n_distinct(t_nodes$loc_num))),
    perc_cover = round(n_sites_cover / n_sites * 100, 0)
  ) %>% filter(focal != coveredby)
  
  comp_df
  
  ggplot() +
    geom_raster(data=comp_df, aes(x=focal, y=coveredby, fill=perc_cover)) +
    scale_fill_gradient(limits = c(0,100)) +
    geom_text(data=comp_df, 
              aes(x=focal, y=coveredby, label=perc_cover), color = "red", size=5) +
    # scale_fill_gradient(high = "#132B43", low = "#56B1F7", limits = c(0,100)) +
    # scale_fill_distiller(palette ="Blues", direction = 1) +
    theme_bw() + theme(panel.grid.major = element_blank()) +
    xlab("Focal") + ylab("Covered by") + labs(fill="% coverage")
  
  # SAVE ##
  ggsave(paste0("figures/site_coverage_compare_", season, "_iba_hex_10km.png"), width=5.5, height = 4)
  
  ## --------------------------------------------------------------------------
  ## Evaluate coverage within latitudinal bands -------------------------------
  ## --------------------------------------------------------------------------
  
  ## Calculate differences in metric ranks btwn datasets
  m_df <- as.data.frame(m_nodes) %>% mutate(datatype="metal")
  c_df <- as.data.frame(c_nodes) %>% mutate(datatype="color")
  t_df <- as.data.frame(t_nodes) %>% mutate(datatype="tracking")
  
  allnodes <- bind_rows(m_df, c_df, t_df) %>% 
    dplyr::select(loc_num, SitRecID, datatype, site_poly, geometry) %>% 
    mutate(longitude = unlist(purrr::map(geometry, 1)),
           latitude = unlist(purrr::map(geometry, 2)))
  
  ## split sites into equal-sized/length latitudinal groups
  # allnodes$lat_group <- as.numeric(cut_number(allnodes$latitude, 4)) # equal-size
  allnodes$lat_group <- as.numeric(cut(allnodes$latitude, 4)) # equal-length
  
  ggplot() +
    geom_histogram(data=allnodes, aes(x=latitude, fill=lat_group, group=lat_group))
  
  ### classify which datatypes 'discovered' each site
  
  site_coince <- allnodes %>% 
    group_by(loc_num, lat_group) %>% 
    summarise(
      n_discover = n_distinct(datatype),
      datatypes  = paste(sort(unique(datatype)), collapse = " "),
      metal      = str_detect(datatypes, pattern = "metal"),
      color      = str_detect(datatypes, pattern = "color"),
      tracking   = str_detect(datatypes, pattern = "tracking")
  ) %>% mutate(
    cat_discover = case_when( ## which sites 'discovered' by each datatype
      metal == TRUE    & (color == FALSE & tracking == FALSE) ~ "M",
      color == TRUE    & (metal == FALSE & tracking == FALSE) ~ "C",
      tracking == TRUE & (metal == FALSE & color == FALSE) ~ "T",
      (metal == TRUE & color == TRUE)    & tracking == FALSE ~ "MC",
      (metal == TRUE & tracking == TRUE) & color == FALSE ~ "MT",
      (color == TRUE & tracking == TRUE) & metal == FALSE ~ "CT",
      TRUE ~ "MCT"
    )
  ) %>% mutate(
    # cat_discover = factor(cat_discover, levels = c("M", "C" , "T", "MC", "MT", "CT", "MCT"))
    cat_discover = factor(cat_discover, levels = c("M", "MC", "C" ,"CT", "T", "MT", "MCT")),
  ) %>% 
    dplyr::select(-datatypes, -metal, -color, -tracking)
  
  ## add geometry back in 
  site_coince %<>% right_join(allnodes[, c("loc_num", "geometry")])
  
  ## save combined site/node layer
  
  saveRDS(site_coince, "data/analysis/networks/allnodes_all_iba_hex_10km.rds.rds")
  
  
  ## summarize by latitude grouping
  latgrp_tot <- site_coince %>% group_by(lat_group) %>% 
    summarise(
      n_sites_tot = n_distinct(loc_num),
      )
  
  latgrp_summ <- site_coince %>% group_by(lat_group, cat_discover) %>% 
    summarise(
      n_sites_cat = n_distinct(loc_num),
    ) %>% left_join(latgrp_tot) %>% 
    mutate(
      perc_sites = (n_sites_cat / n_sites_tot) * 100,
      lat_group = factor(lat_group, levels = c("1", "2" , "3", "4")),
    )
  
  ggplot() +
    geom_col(data=latgrp_summ, aes(x=perc_sites, y=lat_group, fill=cat_discover)) +
    scale_fill_manual(
      values=c("M"="red", "MC"="orange", "C"="yellow", "CT"="green", "T"= "blue", 
               "MT"="purple", "MCT"="black"),
      labels=c("M"="Metal", "MC"="Metal/Color", "C"="Color", "CT"="Color/Track", "T"= "Track","MT"="Metal/Track", "MCT"="All")
    ) +
    theme_bw() + 
    theme(
      panel.grid = element_blank()
    )
  
  ggsave("figures/latgroup_sitediscover_perc.png", width=5, height=7)

  
  ## maps  --------------------------------------------------------------------
  
  mapit <- site_coince
  
  ## optionally remove 'all' category
  # mapit <- mapit %>% filter(cat_discover != "MCT") %>% 
  #   mutate(cat_discover = factor(
  #     cat_discover, 
  #     levels = c("M", "C" , "T", "MC", "MT", "CT")))
  
  # get world map
  wmap <- rworldmap::getMap(resolution="high")
  
  ## Lambert EA: EPSG:3035 ## Lambert CConic: EPSG:3034
  # rams_sf_prj <- st_as_sf(rams) %>% st_transform(crs = "EPSG:3035")
  # iba_sf_prj <- st_as_sf(iba) %>% st_transform(crs = "EPSG:3035")
  wmap_prj <- st_as_sf(wmap) %>% st_transform(crs = "EPSG:3035")
  
  ## Map of all sites -----------------------------------------------------------
  bbox_prj <- st_bbox( ## EAF bbox
    c(xmin = -16, xmax = 21,
      ymin = 8, ymax = 65.2), crs = 4326) %>% 
    st_as_sfc() %>% st_transform("EPSG:3035") %>% st_bbox()
  # bbox_prj <- st_bbox( ## English Channel ---------
  #   c(xmin = -4, xmax = 10,
  #     ymin = 46, ymax = 54), crs = 4326) %>% 
  #   st_as_sfc() %>% st_transform("EPSG:3035") %>% st_bbox()
  # bbox_prj <- st_bbox( ## Iberia ------------------
  #   c(xmin = -10, xmax = 3,
  #     ymin = 35, ymax = 44), crs = 4326) %>% 
  #   st_as_sfc() %>% st_transform("EPSG:3035") %>% st_bbox()
  
  ## save theme to re-use
  maptheme <- 
    theme(
      # plot.background = element_rect(fill = "white"),
      panel.background = element_rect(fill = "white"),
      # axis.text = element_blank(),
      # axis.ticks = element_blank(),
      # panel.grid.major = element_line(colour="grey85"),
      panel.grid.major = element_line(colour="grey60", size=.2),
      # panel.grid.major = element_blank(),# panel.grid.minor = element_blank(),
      panel.border = element_rect(colour = "black", fill=NA, size=.5),
      legend.key=element_blank()
    )
  
  ### 
  node_prj <- st_as_sf(mapit) %>% st_transform("EPSG:3035")
  
  map <- ggplot() +
    geom_sf(data = wmap_prj, fill = "grey70", color = "white", size=0.2) +
    geom_sf(data = node_prj,  
            aes(color = cat_discover)) +
    coord_sf(xlim = c(bbox_prj[1], bbox_prj[3]), 
             ylim = c(bbox_prj[2], bbox_prj[4]), 
             expand = T) +
    scale_size(range = c(0.01, 2)) +
    scale_color_manual(
      values=c("M"="red", "MC"="orange", "C"="yellow", "CT"="green", "T"= "blue", 
               "MT"="purple", "MCT"="black"),
      labels=c("M"="Metal", "MC"="Metal/Color", "C"="Color", "CT"="Color/Track", "T"= "Track","MT"="Metal/Track", "MCT"="All")
    ) +
    maptheme + 
    labs(color="Data type")
  
  map
  
  ## Save ##
  ggsave(paste0("figures/networks/site_discover_compare_", season, "_iba_hex_10kmX.png"), 
         plot=map, width=5, height = 6)
}