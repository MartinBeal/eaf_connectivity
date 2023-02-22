## Compare rank differences between networks of differing data types

pacman::p_load(sfnetworks, sf, dplyr, magrittr, mapview, ggplot2, patchwork)


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
  
  # write.csv(
  #   comp_df, paste0(
  #     "data/analysis/summaries/compare_networks_pairwise_deg10_node_coverage_", season, ".csv",
  #     row.names = F)
  
   ## ----------------------------------------------------------------------------
  ## Combine all three ----------------------------------------------------------
  ## ----------------------------------------------------------------------------
  
  ## Calculate differences in metric ranks btwn datasets
  m_df <- as.data.frame(m_nodes)
  c_df <- as.data.frame(c_nodes)
  t_df <- as.data.frame(t_nodes)
  
  
  ## which sites are in top 10% for degree
  m_df %<>% dplyr::select(
    loc_num, n_id, n_obs, degree, degree_rank, between, btwn_rank, geometry) %>% 
    rename(n_id.m = n_id, n_obs.m = n_obs, degree.m = degree, 
           degree_rank.m = degree_rank, between.m = between, 
           btwn_rank.m = btwn_rank) %>% 
    mutate(
      m_deg10  = ifelse(degree.m >= quantile(degree.m, .9), TRUE, FALSE),
      m_btwn10 = ifelse(between.m >= quantile(between.m, .9), TRUE, FALSE)
    )
  
  c_df %<>% dplyr::select(
    loc_num, n_id, n_obs, degree, degree_rank, between, btwn_rank, geometry) %>% 
    rename(n_id.c = n_id, n_obs.c = n_obs, degree.c = degree, 
           degree_rank.c = degree_rank, between.c = between, 
           btwn_rank.c = btwn_rank) %>% 
    mutate(
      c_deg10  = ifelse(degree.c >= quantile(degree.c, .9), TRUE, FALSE),
      c_btwn10 = ifelse(between.c >= quantile(between.c, .9), TRUE, FALSE)
    )
  
  t_df %<>% dplyr::select(
    loc_num, n_id, n_obs, degree, degree_rank, between, btwn_rank, geometry) %>% 
    rename(n_id.t = n_id, n_obs.t = n_obs, degree.t = degree, 
           degree_rank.t = degree_rank, between.t = between, 
           btwn_rank.t = btwn_rank) %>% 
    mutate(
      t_deg10  = ifelse(degree.t >= quantile(degree.t, .9), TRUE, FALSE),
      t_btwn10 = ifelse(between.t >= quantile(between.t, .9), TRUE, FALSE)
    )
  
  ## merge 
  onetwothree <- plyr::join_all(
    list(m_df, c_df, t_df), 
    by='loc_num', type='full'
  )
  
  onetwothree %<>% mutate(
    m_deg10 = ifelse(is.na(m_deg10), FALSE, m_deg10),
    c_deg10 = ifelse(is.na(c_deg10), FALSE, c_deg10),
    t_deg10 = ifelse(is.na(t_deg10), FALSE, t_deg10),
    m_btwn10 = ifelse(is.na(m_btwn10), FALSE, m_btwn10),
    c_btwn10= ifelse(is.na(c_btwn10), FALSE, c_btwn10),
    t_btwn10= ifelse(is.na(t_btwn10), FALSE, t_btwn10)
  )
  
  ## how similar are ranks btwn metrics for each data type ----------------------
  mplot <- ggplot() + 
    geom_point(data = onetwothree, aes(degree_rank.m, btwn_rank.m)) + theme_bw()+
    ggtitle("Metal") + xlab("Degree rank") + ylab("Betweenness rank")
  cor.test(onetwothree$degree_rank.m, onetwothree$btwn_rank.m)
  
  cplot <- ggplot() + 
    geom_point(data = onetwothree, aes(degree_rank.c, btwn_rank.c)) + theme_bw()  +
    ggtitle("Color") + xlab("Degree rank") + ylab("Betweenness rank")
  cor.test(onetwothree$degree_rank.c, onetwothree$btwn_rank.c) 
  
  tplot <- ggplot() + 
    geom_point(data = onetwothree, aes(degree_rank.t, btwn_rank.t)) + theme_bw() +
    ggtitle("Tracking") + xlab("Degree rank") + ylab("Betweenness rank")
  cor.test(onetwothree$degree_rank.t, onetwothree$btwn_rank.t)
  
  combplot <- mplot / cplot / tplot
  
  ggsave(paste0("figures/deg_btwn_compare_", season,".png"), width = 3, height=8)
  
  
  ## how many datatypes 'discovered' use of each site? --------------------------
  
  ## only sites in top 10%
  onetwothree %<>% group_by(site_poly) %>% mutate(
    n_discover = sum(!is.na(n_id.m), !is.na(n_id.c), !is.na(n_id.t)),
    cat_discover = case_when( ## which sites 'discovered' by each datatype
      !is.na(n_id.m) & (is.na(n_id.c) & is.na(n_id.t)) ~ "M",
      !is.na(n_id.c) & (is.na(n_id.m) & is.na(n_id.t)) ~ "C",
      !is.na(n_id.t) & (is.na(n_id.m) & is.na(n_id.c)) ~ "T",
      (!is.na(n_id.m) & !is.na(n_id.c)) & is.na(n_id.t) ~ "MC",
      (!is.na(n_id.m) & !is.na(n_id.t)) & is.na(n_id.c) ~ "MT",
      (!is.na(n_id.c) & !is.na(n_id.t)) & is.na(n_id.m) ~ "CT",
      TRUE ~ "MCT"
    ),
    cat_deg10 = case_when( ## only sites in top 10%
      m_deg10  == T & (c_deg10 == F & t_deg10 == F) ~ "M",
      c_deg10  == T & (m_deg10 == F & t_deg10 == F) ~ "C",
      t_deg10  == T & (m_deg10 == F & c_deg10 == F) ~ "T",
      (m_deg10 == T & c_deg10 == T) & t_deg10 == F ~ "MC",
      (m_deg10 == T & t_deg10 == T) & c_deg10  == F ~ "MT",
      (c_deg10 == T & t_deg10 == T) & m_deg10 == F ~ "CT",
      (c_deg10 == T & t_deg10 == T & m_deg10 == T) ~ "MCT",
      TRUE ~ "NA" # sites not in 10% get NA
    ),
    cat_btwn10 = case_when( ## only sites in top 10%
      m_btwn10  == T & (c_btwn10 == F & t_btwn10 == F) ~ "M",
      c_btwn10  == T & (m_btwn10 == F & t_btwn10== F) ~ "C",
      t_btwn10  == T & (m_btwn10 == F & c_btwn10 == F) ~ "T",
      (m_btwn10 == T & c_btwn10 == T) & t_btwn10 == F ~ "MC",
      (m_btwn10 == T & t_btwn10 == T) & c_btwn10  == F ~ "MT",
      (c_btwn10 == T & t_btwn10 == T) & m_btwn10 == F ~ "CT",
      (c_btwn10 == T & t_btwn10 == T & m_btwn10 == T) ~ "MCT",
      TRUE ~ "NA" # sites not in 10% get NA
    ),
    cat_discover = factor(cat_discover, levels = c("M", "C" , "T", "MC", "MT", "CT", "MCT")),
    cat_deg10 = factor(cat_deg10, levels = c("M", "C" , "T", "MC", "MT", "CT", "MCT")),
    cat_btwn10 = factor(cat_btwn10, levels = c("M", "C" , "T", "MC", "MT", "CT", "MCT"))
  )
  
  # st_as_sf(onetwothree) %>% mapview(zcol="n_discover")
  # st_as_sf(onetwothree) %>% mapview(zcol="cat_discover")
  # st_as_sf(onetwothree) %>% mapview(zcol="cat_deg10")
  st_as_sf(onetwothree) %>% filter(!is.na(cat_deg10)) %>% mapview(zcol="cat_deg10")
  st_as_sf(onetwothree) %>% filter(!is.na(cat_btwn10)) %>% mapview(zcol="cat_btwn10")
  
  ## Summarize
  tot_sites     <- n_distinct(onetwothree$site_poly)
  discover_summ <- onetwothree %>% 
    group_by(cat_discover) %>% 
    summarise(
      n_sites = n_distinct(site_poly),
      perc_tot = round(n_sites / tot_sites * 100, 1) # % sites found by each data type combo
    )
  
  ## Summarize how many 
  global_compare$perc_tot_nodes <- round(c(NA,
    sum(discover_summ$n_sites[c(1,4,5,6)])/tot_sites,
    sum(discover_summ$n_sites[c(2,4,6,7)])/tot_sites,
    sum(discover_summ$n_sites[c(3,5,6,7)])/tot_sites
  )*100,1)
  
  ## % of ALL sites only found by each datatype (i.e. by only one type)
  global_compare$perc_tot_unique_nodes <- c(NA, discover_summ$perc_tot[1:3])
  
  ## % of sites discovered by each datatype unique to it (only found by it)
  global_compare$perc_rel_unique_nodes <- round(c(NA,
    discover_summ$n_sites[1]/sum(discover_summ$n_sites[c(1,4,5,6)]),
    discover_summ$n_sites[2]/sum(discover_summ$n_sites[c(2,4,6,7)]),
    discover_summ$n_sites[3]/sum(discover_summ$n_sites[c(3,5,6,7)])
  )*100,1)
  
  ## SAVE ##
  write.csv(
    global_compare, 
    paste0("data/analysis/summaries/compare_networks_", season, "_global_metrics.csv"),
    row.names = F)
  
  
  ## maps  --------------------------------------------------------------------
  
  mapit <- onetwothree
  
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
      plot.background = element_rect(fill = "white"),
      panel.background = element_rect(fill = "white"),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
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
  
  
  ### Zoom in on certain regions ------------------------------------------------
  ## English Channel ---------
  bbox2_prj <- st_bbox(
    c(xmin = -4, xmax = 10,
      ymin = 46, ymax = 54), crs = 4326) %>%
    st_as_sfc() %>% st_transform("EPSG:3035") %>% st_bbox()
  
  map2 <- map +
    coord_sf(xlim = c(bbox2_prj[1], bbox2_prj[3]), 
             ylim = c(bbox2_prj[2], bbox2_prj[4]), 
             expand = T)
  
  ## Iberia ------------------
  bbox3_prj <- st_bbox(
    c(xmin = -10, xmax = 3,
      ymin = 35, ymax = 44), crs = 4326) %>%
    st_as_sfc() %>% st_transform("EPSG:3035") %>% st_bbox()
  
  map3 <- map +
    coord_sf(xlim = c(bbox3_prj[1], bbox3_prj[3]), 
             ylim = c(bbox3_prj[2], bbox3_prj[4]), 
             expand = T)
  
  
  ## Save ##
  ggsave(paste0("figures/networks/site_discover_compare_", season, "_iba10km.png"), 
         plot=map, width=5, height = 6)
  ggsave(paste0("figures/networks/site_discover_compare_", season, "_iba10km_EChannel.png"), plot=map2, width=7, height = 6)
  ggsave(paste0("figures/networks/site_discover_compare_", season, "_iba10km_Iberia.png"), plot=map3, width=7, height = 6)
  
  
  ## ----------------------------------------------------------------------------
  ## Map of 10% sites -----------------------------------------------------------
  ## ----------------------------------------------------------------------------
  
  ## Degree - top 10% 
  mapit2 <- st_as_sf(onetwothree) %>% filter(!is.na(cat_deg10))
  
  node_prj <- st_as_sf(mapit2) %>% st_transform("EPSG:3035")
  
  bbox_prj <- st_bbox( ## EAF bbox
    c(xmin = -16, xmax = 21,
      ymin = 8, ymax = 65.2), crs = 4326) %>%
    st_as_sfc() %>% st_transform("EPSG:3035") %>% st_bbox()
  
  
  ## map degree -----------------------------------------------------------------
  map_d <- ggplot() +
    geom_sf(data = wmap_prj, fill = "grey70", color = "white", size=0.2) +
    geom_sf(data = node_prj,  
            aes(color = cat_deg10)) +
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
  
  ### Zoom in on certain regions ------------------------------------------------
  
  ## English Channel ---------
  bbox2_prj <- st_bbox( ## English Channel
    c(xmin = -4, xmax = 10,
      ymin = 46, ymax = 54), crs = 4326) %>%
    st_as_sfc() %>% st_transform("EPSG:3035") %>% st_bbox()
  
  map_d2 <- map_d +
    coord_sf(xlim = c(bbox2_prj[1], bbox2_prj[3]), 
             ylim = c(bbox2_prj[2], bbox2_prj[4]), 
             expand = T)
  
  ## Iberia ---------
  bbox3_prj <- st_bbox( ## Iberia
    c(xmin = -10, xmax = 3,
      ymin = 35, ymax = 44), crs = 4326) %>%
    st_as_sfc() %>% st_transform("EPSG:3035") %>% st_bbox()
  
  map_d3 <- map_d +
    coord_sf(xlim = c(bbox3_prj[1], bbox3_prj[3]), 
             ylim = c(bbox3_prj[2], bbox3_prj[4]), 
             expand = T)
  
  
  
  ## Save ##
  ggsave(paste0("figures/networks/deg10_datatype_compare_", season, "_iba10km.png"), plot=map_d, width=5, height = 6)
  ggsave(paste0("figures/networks/deg10_datatype_compare_", season, "_iba10km_EChannel.png"), plot=map_d2, width=7, height = 6)
  ggsave(paste0("figures/networks/deg10_datatype_compare_", season, "_iba10km_Iberia.png"), plot=map_d3, width=7, height = 6)
  
  
  ## map betweenness ------------------------------------------------------------
  mapit3 <- st_as_sf(onetwothree) %>% filter(!is.na(cat_btwn10))
  node_prj <- st_as_sf(mapit3) %>% st_transform("EPSG:3035")
  
  map_b <- ggplot() +
    geom_sf(data = wmap_prj, fill = "grey70", color = "white", size=0.2) +
    geom_sf(data = node_prj,  
            aes(color = cat_btwn10)) +
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
  
  ### Zoom in on certain regions ------------------------------------------------
  
  ## English Channel ---------
  bbox2_prj <- st_bbox( ## English Channel
    c(xmin = -4, xmax = 10,
      ymin = 46, ymax = 54), crs = 4326) %>%
    st_as_sfc() %>% st_transform("EPSG:3035") %>% st_bbox()
  
  map_b2 <- map_b +
    coord_sf(xlim = c(bbox2_prj[1], bbox2_prj[3]), 
             ylim = c(bbox2_prj[2], bbox2_prj[4]), 
             expand = T)
  
  ## Iberia ---------
  bbox3_prj <- st_bbox( ## Iberia
    c(xmin = -10, xmax = 3,
      ymin = 35, ymax = 44), crs = 4326) %>%
    st_as_sfc() %>% st_transform("EPSG:3035") %>% st_bbox()
  
  map_b3 <- map_b +
    coord_sf(xlim = c(bbox3_prj[1], bbox3_prj[3]), 
             ylim = c(bbox3_prj[2], bbox3_prj[4]), 
             expand = T)
  
  
  ## Save ##
  ggsave(paste0("figures/networks/btwn10_datatype_compare_", season, "_iba10km.png"), plot=map_b, width=5, height = 6)
  ggsave(paste0("figures/networks/btwn10_datatype_compare_", season, "_iba10km_EChannel.png"), plot=map_b2, width=7, height = 6)
  ggsave(paste0("figures/networks/btwn10_datatype_compare_", season, "_iba10km_Iberia.png"), plot=map_b3, width=7, height = 6)
  
  
  
  ## how many of top 10% sites are shared by btwnness and degree? ---------------
  xx <- filter(onetwothree, !is.na(cat_btwn10) | !is.na(cat_deg10))
  
  ## % agreement
  100 - round(sum(is.na(xx$cat_btwn10)) / nrow(xx)*100, 1) ## % btwn cover deg
  100 - round(sum(is.na(xx$cat_deg10)) / nrow(xx)*100, 1)  ## % deg cover btwn
  
  
  ## DEGREE --- Correspondence matrix of top 10% sites found by data types ------
  comp_df <- expand.grid(
    c("metal", "color", "tracking"),c("metal", "color", "tracking"))
  comp_df %<>% dplyr::select(Var2, Var1) %>% rename(typeA = Var2, typeB=Var1)
  
  ## how well does data type B cover data type A? 
  # eg. how many top sites identified from tracking are id'd from color resights?
  m_deg10 <- subset(m_df, m_df$m_deg10 == T)
  c_deg10 <- subset(c_df, c_df$c_deg10 == T)
  t_deg10 <- subset(t_df, t_df$t_deg10 == T)
  
  comp_df$n_sites_cover <- c(
    n_distinct(m_nodes$site_poly),
    sum(c_deg10$site_poly %in% m_deg10$site_poly),
    sum(t_deg10$site_poly %in% m_deg10$site_poly),
    sum(m_deg10$site_poly %in% c_deg10$site_poly),
    n_distinct(c_nodes$site_poly),
    sum(t_deg10$site_poly %in% c_deg10$site_poly),
    sum(m_deg10$site_poly %in% t_deg10$site_poly),
    sum(c_deg10$site_poly %in% t_deg10$site_poly),
    n_distinct(t_deg10$site_poly)
  )
  
  comp_df %<>% mutate(
    n_sites = ifelse(typeA == "metal", n_distinct(m_deg10$site_poly),
                     ifelse(typeA == "color", n_distinct(c_deg10$site_poly),
                            n_distinct(t_deg10$site_poly))),
    perc_cover = round(n_sites_cover / n_sites *100,0)
  ) %>% filter(typeA != typeB)
  comp_df
  
  
  ## Pairwise grid of datatype, top site coverage
  ggplot() +
    geom_raster(data=comp_df, aes(x=typeA, y=typeB, fill=perc_cover)) +
    scale_fill_gradient(limits = c(0,100)) +
    geom_text(data=comp_df, 
              aes(x=typeA, y=typeB, label=perc_cover), color = "red", size=5) + ## add values to grid
    theme_bw() + theme(panel.grid.major = element_blank()) + ggtitle("Degree") +
    xlab("Covered") + ylab("Covering") + labs(fill="% coverage")
  
  ## SAVE ##
  ggsave(paste0("figures/deg10_coverage_compare_", season, "_iba10km.png"), width=5.5, height = 4)
  
  comp_df <- expand.grid(
    c("metal", "color", "tracking"),c("metal", "color", "tracking"))
  comp_df %<>% dplyr::select(Var2, Var1) %>% rename(typeA = Var2, typeB=Var1)
  
  
  ## BETWEENNESS --- Correspondence matrix of top 10% sites found by data types -
  comp_df2 <- expand.grid(
    c("metal", "color", "tracking"),c("metal", "color", "tracking"))
  comp_df2 %<>% dplyr::select(Var2, Var1) %>% 
    rename(typeA = Var2, typeB=Var1)## how well does data B cover data A? 
  
  # eg. how many sites identified from tracking are id'd from color resights?
  m_btwn10 <- subset(m_df, m_df$m_btwn10 == T)
  c_btwn10 <- subset(c_df, c_df$c_btwn10 == T)
  t_btwn10 <- subset(t_df, t_df$t_btwn10 == T)
  
  comp_df2$n_sites_cover <- c(
    n_distinct(m_btwn10$site_poly),
    sum(c_btwn10$site_poly %in% m_btwn10$site_poly),
    sum(t_btwn10$site_poly %in% m_btwn10$site_poly),
    sum(m_btwn10$site_poly %in% c_btwn10$site_poly),
    n_distinct(c_btwn10$site_poly),
    sum(t_btwn10$site_poly %in% c_btwn10$site_poly),
    sum(m_btwn10$site_poly %in% t_btwn10$site_poly),
    sum(c_btwn10$site_poly %in% t_btwn10$site_poly),
    n_distinct(t_btwn10$site_poly)
  )
  
  comp_df2 %<>% mutate(
    n_sites = ifelse(typeA == "metal", n_distinct(m_btwn10$site_poly),
                     ifelse(typeA == "color", n_distinct(c_btwn10$site_poly),
                            n_distinct(t_btwn10$site_poly))),
    perc_cover = round(n_sites_cover / n_sites *100,0)
  ) %>% filter(typeA != typeB)
  comp_df2
  
  ggplot() +
    geom_raster(data=comp_df2, aes(x=typeA, y=typeB, fill=perc_cover)) +
    scale_fill_gradient(limits = c(0,100)) +
    geom_text(data=comp_df2, 
              aes(x=typeA, y=typeB, label=perc_cover), color = "red", size=5) +
    theme_bw() + theme(panel.grid.major = element_blank()) + ggtitle("Betweenness") +
    xlab("Covered") + ylab("Covering") + labs(fill="% coverage")
  
  ## SAVE ##
  ggsave(paste0("figures/btwn10_coverage_compare_", season, "_iba10kmX.png"), width=5.5, height = 4)
  

}
