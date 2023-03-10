## Compare rank differences between networks of differing data types

pacman::p_load(sfnetworks, sf, dplyr, magrittr, mapview, ggplot2, patchwork)

m_net <- readRDS("data/analysis/networks/metal_all_iba10km_poly.rds")
c_net <- readRDS("data/analysis/networks/color_all_iba10km_poly.rds")
t_net <- readRDS("data/analysis/networks/trax_all_iba10km_poly.rds")

## interactive map it
anet <- m_net
nodesf <- anet %>% activate("nodes") %>% sf::st_as_sf()
edgesf <- anet %>% activate("edges") %>% sf::st_as_sf()
mapview::mapview(nodesf, zcol="n_id")
mapview::mapview(edgesf)

m_nodes <- st_as_sf(m_net, "nodes")
c_nodes <- st_as_sf(c_net, "nodes")
t_nodes <- st_as_sf(t_net, "nodes")

m_edges <- st_as_sf(m_net, "edges")
c_edges <- st_as_sf(c_net, "edges")
t_edges <- st_as_sf(t_net, "edges")

global_compare <- data.frame(
  data_type = c("metal", "color", "tracking"),
  n_nodes = c(nrow(m_nodes), nrow(c_nodes), nrow(t_nodes)),
  n_edges = c(nrow(m_edges), nrow(c_edges), nrow(t_edges))
)
global_compare


## Correspondence matrix of sites used by data types --------------------------
comp_df <- expand.grid(
  c("metal", "color", "tracking"),c("metal", "color", "tracking"))
comp_df %<>% dplyr::select(Var2, Var1) %>% rename(typeA = Var2, typeB=Var1)

## how well does data B cover data A? 
# eg. how many sites identified from tracking are id'd from color resights?
comp_df$n_sites_cover <- c(
  n_distinct(m_nodes$site_poly),
  sum(c_nodes$site_poly %in% m_nodes$site_poly),
  sum(t_nodes$site_poly %in% m_nodes$site_poly),
  sum(m_nodes$site_poly %in% c_nodes$site_poly),
  n_distinct(c_nodes$site_poly),
  sum(t_nodes$site_poly %in% c_nodes$site_poly),
  sum(m_nodes$site_poly %in% t_nodes$site_poly),
  sum(c_nodes$site_poly %in% t_nodes$site_poly),
  n_distinct(t_nodes$site_poly)
)

comp_df %<>% mutate(
  n_sites = ifelse(typeA == "metal", n_distinct(m_nodes$site_poly),
                   ifelse(typeA == "color", n_distinct(c_nodes$site_poly),
                          n_distinct(t_nodes$site_poly))),
  perc_cover = round(n_sites_cover / n_sites *100,0)
) %>% filter(typeA != typeB)

comp_df

ggplot() +
  geom_raster(data=comp_df, aes(x=typeA, y=typeB, fill=perc_cover)) +
  scale_fill_gradient(limits = c(0,100)) +
  # scale_fill_gradient(high = "#132B43", low = "#56B1F7", limits = c(0,100)) +
  # scale_fill_distiller(palette ="Blues", direction = 1) +
  theme_bw() + theme(panel.grid.major = element_blank()) +
  xlab("Covered") + ylab("Covering") + labs(fill="% coverage")

## SAVE ##
ggsave(paste0("figures/site_coverage_compare_iba10kmX.png"), width=5.5, height = 4)

write.csv(
  comp_df, "data/analysis/summaries/compare_networks_pairwise_deg10_node_coverage.csv", row.names = F)


## Calculate differences in metric ranks btwn datasets ------------------------
m_df <- as.data.frame(m_nodes)
c_df <- as.data.frame(c_nodes)
t_df <- as.data.frame(t_nodes)

## compare two node sets (eg. trax and color)
t_c_df <- full_join(t_df, c_df, by=c("site_poly"))

plot(t_c_df$btwn_rank.x, t_c_df$btwn_rank.y)
plot(t_c_df$degree_rank.x, t_c_df$degree_rank.y)

cor.test(x=t_c_df$btwn_rank.x, y=t_c_df$btwn_rank.y, method = 'spearman')
cor.test(x=t_c_df$degree_rank.x, y=t_c_df$degree_rank.y, method = 'spearman')


## calculate rank differences 
t_c_df <- t_c_df %>% 
  mutate(
    # btwn_rdiff = abs(btwn_rank.x - btwn_rank.y),
    btwn_rdiff = btwn_rank.x - btwn_rank.y,
    btwn_rsum  = sum(btwn_rank.x, btwn_rank.y),
    # degree_rdiff = abs(degree_rank.x - degree_rank.y),
    degree_rdiff = degree_rank.x - degree_rank.y, ## neg. rdiff = more important from trax
    degree_rsum  = sum(degree_rank.x, degree_rank.y)
  )

st_as_sf(t_c_df) %>% mapview::mapview(zcol="btwn_rdiff")
st_as_sf(t_c_df) %>% mapview::mapview(zcol="degree_rdiff")



## Compare only shared sites --------------------------------------------------

pacman::p_load(sfnetworks, sf, dplyr)

m_net <- readRDS("data/analysis/networks/metal_all_iba10km_poly.rds")
c_net <- readRDS("data/analysis/networks/color_all_iba10km_poly.rds")
t_net <- readRDS("data/analysis/networks/trax_all_iba10km_poly.rds")

## get nodes only
m_nodes <- st_as_sf(m_net, "nodes")
c_nodes <- st_as_sf(c_net, "nodes")
t_nodes <- st_as_sf(t_net, "nodes")

## Calculate differences in metric ranks btwn datasets
m_df <- as.data.frame(m_nodes)
c_df <- as.data.frame(c_nodes)
t_df <- as.data.frame(t_nodes)

## three pairwise comparisons to make
onetwo   <- full_join(m_df, c_df, by=c("site_poly")) %>% 
  dplyr::select(site_poly, n_id.x, n_id.y, n_obs.x, n_obs.y, degree.x, degree.y, 
                degree_rank.x, degree_rank.y, between.x, between.y, 
                btwn_rank.x, btwn_rank.y, geometry.x)
onethree <- full_join(m_df, t_df, by=c("site_poly"))%>% 
  dplyr::select(site_poly, n_id.x, n_id.y, n_obs.x, n_obs.y, degree.x, degree.y, 
                degree_rank.x, degree_rank.y, between.x, between.y, 
                btwn_rank.x, btwn_rank.y, geometry.x)
twothree <- full_join(c_df, t_df, by=c("site_poly"))%>% 
  dplyr::select(site_poly, n_id.x, n_id.y, n_obs.x, n_obs.y, degree.x, degree.y, 
                degree_rank.x, degree_rank.y, between.x, between.y, 
                btwn_rank.x, btwn_rank.y, geometry.x)

## Get only sites identified by both datatypes
onetwo$inboth   <- !(is.na(onetwo$degree_rank.x)   | is.na(onetwo$degree_rank.y))
onethree$inboth <- !(is.na(onethree$degree_rank.x) | is.na(onethree$degree_rank.y))
twothree$inboth <- !(is.na(twothree$degree_rank.x) | is.na(twothree$degree_rank.y))

onetwo_f   <- filter(onetwo, inboth == TRUE)
onethree_f <- filter(onethree, inboth == TRUE)
twothree_f <- filter(twothree, inboth == TRUE)

## Metal vs. CR
onetwo_f %<>%
  mutate(
    degree_rank2.x = rank(desc(degree.x), ties.method = "min"),
    degree_rank2.y = rank(desc(degree.y), ties.method = "min"),
    degree_rdiff = abs(degree_rank2.x - degree_rank2.y),
    degree_rdiff2 = degree_rank2.x - degree_rank2.y, ## neg. rdiff = more important from metal
    btwn_rank2.x = rank(desc(between.x), ties.method = "min"),
    btwn_rank2.y = rank(desc(between.y), ties.method = "min"),
    btwn_rdiff = abs(btwn_rank2.x - btwn_rank2.y),
    btwn_rdiff2 = btwn_rank2.x - btwn_rank2.y, ## neg. rdiff = more important from metal
  )

plot(onetwo_f$degree_rank2.x, onetwo_f$degree_rank2.y, pch=19)
plot(onetwo_f$degree.x, onetwo_f$degree.y, pch=19)

st_as_sf(onetwo_f) %>% mapview::mapview(zcol="degree_rdiff2")
st_as_sf(onetwo_f) %>% mapview::mapview(zcol="degree_rdiff")

## Metal vs. trax
onethree_f %<>%
  mutate(
    degree_rank2.x = rank(desc(degree.x), ties.method = "min"),
    degree_rank2.y = rank(desc(degree.y), ties.method = "min"),
    degree_rdiff = abs(degree_rank2.x - degree_rank2.y),
    degree_rdiff2 = degree_rank2.x - degree_rank2.y, ## neg. rdiff = more important from metal
    btwn_rank2.x = rank(desc(between.x), ties.method = "min"),
    btwn_rank2.y = rank(desc(between.y), ties.method = "min"),
    btwn_rdiff = abs(btwn_rank2.x - btwn_rank2.y),
    btwn_rdiff2 = btwn_rank2.x - btwn_rank2.y, ## neg. rdiff = more important from metal
  )

plot(onethree_f$degree_rank2.x, onethree_f$degree_rank2.y, pch=19)
plot(onethree_f$degree.x, onethree_f$degree.y, pch=19)

st_as_sf(onethree_f) %>% mapview::mapview(zcol="degree_rdiff2")
st_as_sf(onethree_f) %>% mapview::mapview(zcol="degree_rdiff")

## CR vs. trax ----------------------------------------------------------------
twothree_f %<>%
  mutate(
    degree_rank2.x = rank(desc(degree.x), ties.method = "min"),
    degree_rank2.y = rank(desc(degree.y), ties.method = "min"),
    degree_rdiff = abs(degree_rank2.x - degree_rank2.y),
    degree_rdiff2 = degree_rank2.x - degree_rank2.y, ## neg. rdiff = more important from metal
    btwn_rank2.x = rank(desc(between.x), ties.method = "min"),
    btwn_rank2.y = rank(desc(between.y), ties.method = "min"),
    btwn_rdiff = abs(btwn_rank2.x - btwn_rank2.y),
    btwn_rdiff2 = btwn_rank2.x - btwn_rank2.y, ## neg. rdiff = more important from metal
  )

plot(twothree_f$degree_rank.x, twothree_f$degree_rank.y, pch=19)
plot(twothree_f$degree_rank2.x, twothree_f$degree_rank2.y, pch=19)
plot(twothree_f$degree.x, twothree_f$degree.y, pch=19)

st_as_sf(twothree_f) %>% mapview::mapview(zcol="degree_rdiff2")
st_as_sf(twothree_f) %>% mapview::mapview(zcol="degree_rdiff")


## ----------------------------------------------------------------------------
## Combine all three ----------------------------------------------------------
## ----------------------------------------------------------------------------

## which sites are in top 10% for degree
m_df %<>% dplyr::select(
  site_poly, n_id, n_obs, degree, degree_rank, between, btwn_rank, geometry) %>% 
  rename(n_id.m = n_id, n_obs.m = n_obs, degree.m = degree, 
         degree_rank.m = degree_rank, between.m = between, 
         btwn_rank.m = btwn_rank) %>% 
  mutate(
    m_deg10  = ifelse(degree.m >= quantile(degree.m, .9), TRUE, FALSE),
    m_btwn10 = ifelse(between.m >= quantile(between.m, .9), TRUE, FALSE)
    )

c_df %<>% dplyr::select(
  site_poly, n_id, n_obs, degree, degree_rank, between, btwn_rank, geometry) %>% 
  rename(n_id.c = n_id, n_obs.c = n_obs, degree.c = degree, 
         degree_rank.c = degree_rank, between.c = between, 
         btwn_rank.c = btwn_rank) %>% 
  mutate(
    c_deg10  = ifelse(degree.c >= quantile(degree.c, .9), TRUE, FALSE),
    c_btwn10 = ifelse(between.c >= quantile(between.c, .9), TRUE, FALSE)
    )

t_df %<>% dplyr::select(
  site_poly, n_id, n_obs, degree, degree_rank, between, btwn_rank, geometry) %>% 
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
  by='site_poly', type='full'
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

ggsave("figures/deg_btwn_compare.png", width = 3, height=8)

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
tot_sites <- n_distinct(onetwothree$site_poly)
discover_summ <- onetwothree %>% group_by(cat_discover) %>% 
  summarise(
    n_sites = n_distinct(site_poly),
    perc_tot = round(n_sites / tot_sites * 100, 1)
  )

global_compare$perc_tot_nodes <- round(c(
  sum(discover_summ$n_sites[c(1,4,5,6)])/tot_sites,
  sum(discover_summ$n_sites[c(2,4,6,7)])/tot_sites,
  sum(discover_summ$n_sites[c(3,5,6,7)])/tot_sites
)*100,1)
global_compare$perc_tot_unique_nodes <- c(discover_summ$perc_tot[1:3])
global_compare$perc_rel_unique_nodes <- round(c(
  discover_summ$n_sites[1]/sum(discover_summ$n_sites[c(1,4,5,6)]),
  discover_summ$n_sites[2]/sum(discover_summ$n_sites[c(2,4,6,7)]),
  discover_summ$n_sites[3]/sum(discover_summ$n_sites[c(3,5,6,7)])
)*100,1)

## SAVE ##
write.csv(global_compare, "data/analysis/summaries/compare_networks_global_metrics.csv", row.names = F)

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
# bbox_prj <- st_bbox( ## English Channel
#   c(xmin = -4, xmax = 10,
#     ymin = 46, ymax = 54), crs = 4326) %>% 
#   st_as_sfc() %>% st_transform("EPSG:3035") %>% st_bbox()
# bbox_prj <- st_bbox( ## Iberia
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

## Save ##
ggsave(paste0("figures/networks/site_discover_compare_iba10kmX.png"), plot=map, width=5, height = 6)
# ggsave(paste0("figures/networks/site_discover_compare_iba10km_zoom2.png"), plot=map, width=7, height = 6)


## Map of 10% sites -----------------------------------------------------------

mapit2 <- st_as_sf(onetwothree) %>% filter(!is.na(cat_deg10))

node_prj <- st_as_sf(mapit2) %>% st_transform("EPSG:3035")

# bbox_prj <- st_bbox( ## EAF bbox
#   c(xmin = -16, xmax = 21,
#     ymin = 8, ymax = 65.2), crs = 4326) %>% 
#   st_as_sfc() %>% st_transform("EPSG:3035") %>% st_bbox()
bbox_prj <- st_bbox( ## English Channel
  c(xmin = -4, xmax = 10,
    ymin = 46, ymax = 54), crs = 4326) %>%
  st_as_sfc() %>% st_transform("EPSG:3035") %>% st_bbox()
# bbox_prj <- st_bbox( ## Iberia
#   c(xmin = -10, xmax = 3,
#     ymin = 35, ymax = 44), crs = 4326) %>% 
#   st_as_sfc() %>% st_transform("EPSG:3035") %>% st_bbox()


## map degree -----------------------------------------------------------------
map2 <- ggplot() +
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

## Save ##
# ggsave(paste0("figures/networks/deg10_datatype_compare_iba10kmX.png"), plot=map2, width=5, height = 6)
ggsave(paste0("figures/networks/deg10_datatype_compare_iba10km_zoom1.png"), plot=map2, width=7, height = 6)
# ggsave(paste0("figures/networks/deg10_datatype_compare_iba10km_zoom2.png"), plot=map2, width=7, height = 6)


## map betweenness ------------------------------------------------------------
mapit3 <- st_as_sf(onetwothree) %>% filter(!is.na(cat_btwn10))
node_prj <- st_as_sf(mapit3) %>% st_transform("EPSG:3035")

map3 <- ggplot() +
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

## Save ##
# ggsave(paste0("figures/networks/btwn10_datatype_compare_iba10kmX.png"), plot=map3, width=5, height = 6)
ggsave(paste0("figures/networks/btwn10_datatype_compare_iba10km_zoom1.png"), plot=map3, width=7, height = 6)
# ggsave(paste0("figures/networks/btwn10_datatype_compare_iba10km_zoom2.png"), plot=map3, width=7, height = 6)

## how many of top 10% sites are shared by btwnness and degree? ---------------
xx <- filter(onetwothree, !is.na(cat_btwn10) | !is.na(cat_deg10))

## % agreement
100 - round(sum(is.na(xx$cat_btwn10)) / nrow(xx)*100, 1) ## % btwn cover deg
100 - round(sum(is.na(xx$cat_deg10)) / nrow(xx)*100, 1)  ## % deg cover btwn


## DEGREE --- Correspondence matrix of top 10% sites found by data types ------
comp_df <- expand.grid(
  c("metal", "color", "tracking"),c("metal", "color", "tracking"))
comp_df %<>% dplyr::select(Var2, Var1) %>% rename(typeA = Var2, typeB=Var1)

## how well does data B cover data A? 
# eg. how many sites identified from tracking are id'd from color resights?
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

ggplot() +
  geom_raster(data=comp_df, aes(x=typeA, y=typeB, fill=perc_cover)) +
  scale_fill_gradient(limits = c(0,100)) +
  theme_bw() + theme(panel.grid.major = element_blank()) + ggtitle("Degree") +
  xlab("Covered") + ylab("Covering") + labs(fill="% coverage")

## SAVE ##
ggsave(paste0("figures/deg10_coverage_compare_iba10kmX.png"), width=5.5, height = 4)

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
  theme_bw() + theme(panel.grid.major = element_blank()) + ggtitle("Betweenness") +
  xlab("Covered") + ylab("Covering") + labs(fill="% coverage")

## SAVE ##
ggsave(paste0("figures/btwn10_coverage_compare_iba10kmX.png"), width=5.5, height = 4)



