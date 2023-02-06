### Compare distributions of edge distances between data types ----------------
## NEED TO DO THIS FOR NODES COMMON TO BOTH NETWORKS

pacman::p_load(sfnetworks, dplyr, sf, mapview, units)

### Season -------------------------------------------------
season <- "all"
# season <- "spring"
# season <- "fall"

## load data
# only IBAs as sites
# trxnet  <- readRDS( ## in sites network
#   paste0("data/analysis/networks/trax_", season, "_iba10km_poly.rds"))
# colnet  <- readRDS( ## in sites network
#   paste0("data/analysis/networks/color_", season, "_iba10km_poly.rds"))
# metnet  <- readRDS( ## in sites network
#   paste0("data/analysis/networks/metal_", season, "_iba10km_poly.rds"))

# IBAs and outside sites
trxnet  <- readRDS( ## in sites network
  paste0("data/analysis/networks/trax_", season, "_iba_hex_10km.rds"))
colnet  <- readRDS( ## in sites network
  paste0("data/analysis/networks/color_", season, "_iba_hex_10km.rds"))
metnet  <- readRDS( ## in sites network
  paste0("data/analysis/networks/metal_", season, "_iba_hex_10km.rds"))

## nodes
trxnodes <- trxnet %>% activate("nodes") %>% st_as_sf()
colnodes <- colnet %>% activate("nodes") %>% st_as_sf()
metnodes <- metnet %>% activate("nodes") %>% st_as_sf()

### Retain only nodes/connections shared btwn networks (intersection)
## Tracking
tsharenet <- st_filter(activate(trxnet, "nodes"), activate(colnet, "nodes"))
# tsharenet <- st_filter(activate(tsharenet, "nodes"), activate(metnet, "nodes"))
## Color 
csharenet <- st_filter(activate(colnet, "nodes"), activate(trxnet, "nodes"))
# csharenet <- st_filter(activate(csharenet, "nodes"), activate(metnet, "nodes"))
## Metal
msharenet <- st_filter(activate(metnet, "nodes"), activate(trxnet, "nodes"))
# msharenet <- st_filter(activate(msharenet, "nodes"), activate(colnet, "nodes"))

## Calculate edge lengths/distances
tsharenet %<>% activate("edges") %>% mutate(
  edge_len = set_units(edge_length(), km)
  # edge_displ  = edge_displacement()
)

csharenet %<>% activate("edges") %>% mutate(
  edge_len = set_units(edge_length(), km)
  # edge_displ  = edge_displacement()
)

## convert edges to sf 
trxedges <- tsharenet %>% activate("edges") %>% st_as_sf() %>% 
  mutate(datatype="Track")
coledges <- csharenet %>% activate("edges") %>% st_as_sf() %>% 
  mutate(datatype="Color")

## maplongest 'direct' connection
coledges[which.max(coledges$edge_len),] %>% mapview()
trxedges[which.max(trxedges$edge_len),] %>% mapview()

combedges <- bind_rows(trxedges, coledges) %>% st_drop_geometry() 

edge_summ <- combedges %>% 
  group_by(datatype) %>% 
  summarise(
    n_edges = n(),
    mn_edge_length = mean(edge_len),
    md_edge_length = median(edge_len),
    sd_edge_length = sd(edge_len)
  )

## data highly non-normal
m1 <- lm(combedges$edge_len ~ combedges$datatype)
qqnorm(residuals(m1))
qqline(residuals(m1))

## boxplot
ggplot() + 
  geom_boxplot(
    data = combedges, aes(x=datatype, y=edge_len, fill=datatype),
    show.legend = FALSE) +
  theme_bw()

ggsave("figures/edge_length_datatype_iba_hex_10km.png", width=5, height=6)


## All three data types -------------------------------------------------------

### Retain only nodes/connections shared btwn networks (intersection)
## Tracking
tsharenet2 <- st_filter(activate(tsharenet, "nodes"), activate(metnet, "nodes"))
## Color 
csharenet2 <- st_filter(activate(csharenet, "nodes"), activate(metnet, "nodes"))
## Metal
msharenet2 <- st_filter(activate(msharenet, "nodes"), activate(colnet, "nodes"))

## Calculate edge lengths/distances
tsharenet2 %<>% activate("edges") %>% mutate(
  edge_len = set_units(edge_length(), km)
  # edge_displ  = edge_displacement()
)

csharenet2 %<>% activate("edges") %>% mutate(
  edge_len = set_units(edge_length(), km)
  # edge_displ  = edge_displacement()
)

msharenet2 %<>% activate("edges") %>% mutate(
  edge_len = set_units(edge_length(), km)
  # edge_displ  = edge_displacement()
)

## convert edges to sf 
trxedges2 <- tsharenet %>% activate("edges") %>% st_as_sf() %>% 
  mutate(datatype="Track")
coledges2 <- csharenet %>% activate("edges") %>% st_as_sf() %>% 
  mutate(datatype="Color")
metedges2 <- msharenet2 %>% activate("edges") %>% st_as_sf() %>% 
  mutate(datatype="Metal")

## maplongest 'direct' connection
trxedges2[which.max(trxedges2$edge_len),] %>% mapview()
coledges2[which.max(coledges2$edge_len),] %>% mapview()
metedges2[which.max(metedges2$edge_len),] %>% mapview()

combedges2 <- bind_rows(trxedges2, coledges2, metedges2) %>% st_drop_geometry() 

edge_summ2 <- combedges2 %>% 
  group_by(datatype) %>% 
  summarise(
    n_edges = n(),
    mn_edge_length = mean(edge_len),
    md_edge_length = median(edge_len),
    sd_edge_length = sd(edge_len)
  )
edge_summ2

write.csv(edge_summ2, "data/analysis/summaries/compare_edge_length_alldatatypes_iba_hex_10km.csv", row.names = F)

## data highly non-normal
m2 <- lm(combedges2$edge_len ~ combedges2$datatype)
qqnorm(residuals(m2))
qqline(residuals(m2))

## boxplot
ggplot() + 
  geom_boxplot(
    data = combedges2, aes(x=datatype, y=edge_len, fill=datatype),
    show.legend = FALSE) +
  theme_bw()

ggsave("figures/edge_length_datatype_alldatatypes_iba_hex_10km.png", width=5, height=6)


## compare basic maps --------------------------------------------------------

trxmap <- ggplot() + geom_sf(data = trxedges2, aes(), alpha=0.4) + theme_void() +
  ggtitle("Tracking")
colmap <- ggplot() + geom_sf(data = coledges2, aes(), alpha=0.4) + theme_void() +
  ggtitle("Color")
metmap <- ggplot() + geom_sf(data = metedges2, aes(), alpha=0.4) + theme_void() +
  ggtitle("Metal")

library(patchwork)

combmap <- colmap + metmap + trxmap

ggsave("figures/networks/compare_edges_shrdnodes.png_iba_hex_10km.png", height=5, width=8)
