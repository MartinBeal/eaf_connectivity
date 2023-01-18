## Map networks and metrics --------------------------------------------------

pacman::p_load(sfnetworks, sf, dplyr, magrittr, ggplot2)

## which network to create
# season <- "all"
season <- "spring"
# season <- "fall"

m_net <- readRDS(
  paste0("data/analysis/networks/metal_", season, "_iba10km_poly.rds"))
c_net <- readRDS(
  paste0("data/analysis/networks/color_", season, "_iba10km_poly.rds"))
t_net <- readRDS(
  paste0("data/analysis/networks/trax_", season, "_iba10km_poly.rds"))

## interactive map it
anet <- m_net
nodesf <- anet %>% activate("nodes") %>% sf::st_as_sf()
edgesf <- anet %>% activate("edges") %>% sf::st_as_sf()
# mapview::mapview(nodesf, zcol="n_id")
# mapview::mapview(edgesf)

m_nodes <- st_as_sf(m_net, "nodes")
c_nodes <- st_as_sf(c_net, "nodes")
t_nodes <- st_as_sf(t_net, "nodes")

m_edges <- st_as_sf(m_net, "edges")
c_edges <- st_as_sf(c_net, "edges")
t_edges <- st_as_sf(t_net, "edges")

## Static map: project for prettier map --------------------------------------

## site polygons
# rams <- raster::shapefile("data/geodata/ramsar/EAF_ramsar.shp")
# iba <- raster::shapefile("data/geodata/ibas/EAF_btgo_IBA/eaf_btgo_iba.shp")

dts <- c("metal", "color", "trax")

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

for(i in seq_along(dts)){
  datatype <- dts[i]
  if(datatype == "metal"){    ## metal ring captures, recaptures, recoveries
    nodes <- m_nodes
    edges <- m_edges
  } else if (datatype == "color"){ ## color ring resightings
    nodes <- c_nodes
    edges <- c_edges
  } else if (datatype == "trax"){  ## tracking locations 
    nodes <- t_nodes
    edges <- t_edges
  }
  ## Lambert EA: EPSG:3035 ## Lambert CConic: EPSG:3034
  node_prj <- nodes %>% st_transform("EPSG:3035")
  edge_prj <- edges %>% st_transform("EPSG:3035")
  
  ## EAF bbox
  bbox_prj <- st_bbox(
    c(xmin = -16, xmax = 21,
      ymin = 8, ymax = 65.2), crs = 4326) %>% 
    st_as_sfc() %>% st_transform("EPSG:3035") %>% st_bbox()
  
  ### Map 1: nodes and edges, weighted by N individuals -----------------------
  map1 <- ggplot() +
    geom_sf(data = wmap_prj, fill = "grey70", color = NA) +
    # geom_sf(data = iba_sf_prj,               ## IBAs
    #         aes(), fill = "dark red", color = NA) +
    geom_sf(data = edge_prj,                 ## edges
            aes(size = n_id), col = "black") + #, alpha = 0.65
    geom_sf(data = wmap_prj, fill = NA, color = "white", size=0.2) +
    geom_sf(data = arrange(node_prj, n_id),  ## n individuals
            aes(col = n_id)) +
    coord_sf(xlim = 
               c(bbox_prj[1], bbox_prj[3]), ylim = c(bbox_prj[2], bbox_prj[4]), 
             expand = T) +
    scale_size(range = c(0.01, 2)) +
    maptheme
  # map1
  
  ## SAVE  ##
  ggsave(
    paste0("figures/networks/", datatype,"_", season, "_iba10kmX.png"), 
    plot=map1, width=5, height = 6)
  # ggsave(paste0("figures/networks/", datatype,"_", season, "_iba10km_polyX.png"), plot=map, width=5, height = 6)
  # ggsave(paste0("figures/networks/", datatype,"_", season, "_iba10km_iceX.png"), plot=map, width=5, height = 6)
  
  ## Map 2:3 network metrics --------------------------------------------------
  map2 <- ggplot() +
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
    geom_sf(data = arrange(node_prj, degree),  ## node degree rank
            aes(color = degree_rank, size = degree_rank)) +
    coord_sf(
      xlim = c(bbox_prj[1], bbox_prj[3]), ylim = c(bbox_prj[2], bbox_prj[4]), 
             expand = T) +
    scale_size(range = c(0.01, 2),   # reverse size scale (for ranks)
               trans = 'reverse') +
    maptheme
  # reverse size scale (for ranks)
    # scale_size(range = c(0.01, 2)) 
  map2
  
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
    geom_sf(data = arrange(node_prj, between),  ## node degree rank
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
    # ggsave(paste0("figures/networks/", datatype,"_", season, "_iba10km_poly_betweeness.png"), plot=map2, width=5, height = 6)
    # ggsave(paste0("figures/networks/", datatype,"_", season, "_iba10km_poly_degree.png"), plot=map2, width=5, height = 6)
    ggsave(
      paste0("figures/networks/", datatype,"_", season, "_iba10km_poly_degreernkX.png"), 
      plot=map2, width=5, height = 6)
    ggsave(
      paste0("figures/networks/", datatype,"_", season, "_iba10km_poly_btwnrnkX.png"), 
      plot=map3, width=5, height = 6)
    
}
