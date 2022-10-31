## TRYING TO SET UP DIRECTED NETWORK

pacman::p_load(dplyr, igraph, stringr, tictoc, tidygraph, sfnetworks, ggplot2, 
               sf, mapview, magrittr, lubridate)

# alldat <- readRDS("data/analysis/cr_tag_ice_wash_merge.rds")
# alldat <- readRDS("data/analysis/cr_tag_ice_wash_merge_districts.rds")
alldat <- readRDS("data/analysis/test_comb_euring_cring_districts.rds")

# fxn for splitting string into columns
source("C:/Users/Martim Bill/Documents/R/source_scripts/str2col.R")

## choose data subset to run --------------------------------------------------

# level <- "site"
level <- "district"
# level <- "country"

## which network to create
season <- "all"
season <- "spring"
season <- "fall"

# netdat <- subset(alldat, combo %in% unique(alldat$combo)[1:100])

if(season == "all"){ ## year-round
  netdat <- alldat ## all individuals
} else { ## separate networks for fall and spring migration 
  alldat$season <- ifelse(month(alldat$date) < 7, "spring", "fall")
  netdat <- subset(alldat, season == season)
}

### Choose which spatial scale to work at -------------------------------------

if(level == "site"){
  netdat$locale <- as.factor(netdat$site)
  netdat$loc_num <- as.numeric(netdat$locale)
} else if(level == "district"){ ## district (county)
  netdat$locale <- as.factor(netdat$dist)
  netdat$loc_num <- as.numeric(netdat$locale)
} else if(level == "country"){ ## Country
  ## Fix some bad country data
  
  netdat$country <- ifelse(netdat$country == "Morroco",
                           "Morocco", netdat$country)
  netdat$country <- ifelse(
    netdat$country %in% c("Scotland", "England",
                          "Channel Islands", "Lincolnshire", "Northern Ireland",
                          "Wales"), "United Kingdom", netdat$country)
  netdat$country <- ifelse(netdat$country == "Faroe Island",
                           "Denmark", netdat$country)
  
  netdat$locale <- as.factor(netdat$country)
  netdat$loc_num <- as.numeric(netdat$locale)
}

netdat %<>% subset(!is.na(loc_num)) # remove sightings w/ NAs in site info

# create reference table of sites w/ obs
site_summ <- netdat %>% group_by(loc_num, locale) %>% 
  summarise(latitude = first(latitude), longitude = first(longitude))

### Quickly get country centroids 

library(rgeos)
library(rworldmap)
# get world map
wmap <- getMap(resolution="high")

if(level == "country"){ ## Country
  centroids <- gCentroid(wmap, byid=TRUE)
  
  ccoords <- as.data.frame(centroids)
  
  ccoords$country <- row.names(ccoords)
  ccoords <- rename(ccoords, latitude=y, longitude=x)
  
  site_summ <- left_join(
    site_summ[,-c(3:4)], ccoords, by = c("locale"="country"))
}


## Edge list ------------------------------------------------------------------

### Directed -- only temporally ordered connections (optionally more criteria)

## try vector-shifting and ifelse

netdats <- select(netdat, combo, date, loc_num) %>% 
  mutate(year = year(date))

# netdats <- subset(netdats, combo == unique(netdats$combo)[1])

# from <- netdats[1:nrow(netdats)-1, ] %>%
#   rename(date_from=date, from=loc_num, year_from=year)
# to   <- netdats[-1, ] %>%
#   rename(date_to=date, to=loc_num, year_to=year)
# 
# netdir <- cbind(from, to[, !names(to) %in% c("combo")])

## apply per individual
netdir_list <- lapply(split(netdats, netdats$combo), function(x){
  
  from <- x[1:nrow(x)-1, ] %>% 
    rename(date_from=date, from=loc_num, year_from=year)
  to   <- x[-1, ] %>% 
    rename(date_to=date, to=loc_num, year_to=year)
  netdir <- cbind(from, to[, !names(to) %in% c("combo")])
  return(netdir)
})

netdir <- do.call(rbind, netdir_list)


## filter to only movements within a single year (same migration season)
netdir_f <- subset(netdir, year_from == year_to)

## remove self connections
noself <- netdir_f[-which(netdir_f$from == netdir_f$to), ]

## combine sites into single variable for summarizing
noself$sitecomb <- paste(noself$from, noself$to)

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
nodelist <- netdat %>% group_by(loc_num) %>% 
  summarise(
    n_id = n_distinct(combo),
    n_obs  = n()
  )

nodelist <- nodelist %>% 
  left_join(site_summ) %>% 
  sf::st_as_sf(coords = c("longitude", "latitude"), 
               crs = 4326, agr = "constant")


## filter out weak connections *****GIVE AN ERROR ??*********
# loc_nums <- unique(netdat$loc_num)
# nodelist %<>% filter(n_id > 1)
# out_locs <- loc_nums[!loc_nums %in% unique(nodelist$loc_num)]
# 
# edgelist %<>% filter(
#   !from %in% out_locs & !to %in% out_locs)
# 
# edgelist %<>% filter(n_id > 1)
# nodelist %<>% filter(n_id > 1)

## Convert to sfnetwork

# netsf <- sfnetwork(nodelist, edgelist, node_key = "dist", directed = F)
netsf <- sfnetwork(nodelist, edgelist, node_key = "dist", directed = T, 
                   edges_as_lines = TRUE)

## if it gives an error try running via igraph
# graph_from_data_frame(
#   d = edgelist, 
#   vertices = nodelist,
#   directed = F)

# netsf <- netsf %>%
#   activate("nodes") %>%
#   filter(loc_num %in% unique(c(edgelist$from, edgelist$to)))

## filter out sites (and links) w/ few observations
# netsf <- netsf %>%
#   activate("nodes") %>%
#   filter(n_id > 4)

## interactive map it
# nodesf <- netsf %>% activate("nodes") %>% sf::st_as_sf()
# netsf %>% activate("edges") %>% sf::st_as_sf() %>% mapview::mapview() + 
#   mapview::mapview(nodesf, zcol="n_id")


## Static map: project for prettier map --------------------------------------

## Lambert EA: EPSG:3035 ## Lambert CConic: EPSG:3034
wmap_prj <- st_as_sf(wmap) %>% st_transform(crs = "EPSG:3035")
edge_prj <- st_as_sf(netsf, "edges") %>% st_transform(crs = "EPSG:3035")
node_prj <- st_as_sf(netsf, "nodes") %>% st_transform("EPSG:3035")
bbox_prj <- st_bbox(edge_prj)

### SEEMS TO BE A BUG (plot.background doesn't work when data are projected)
map <- ggplot() +
  geom_sf(data = wmap_prj, fill = "grey80", color = NA) +
  geom_sf(data = edge_prj, 
          aes(size = n_id), col = "black") + #, alpha = 0.65
  geom_sf(data = wmap_prj, fill = NA, color = "white") +
  geom_sf(data = node_prj, 
          aes(col = n_id)) +
  coord_sf(xlim = 
             c(bbox_prj[1], bbox_prj[3]), ylim = c(bbox_prj[2], bbox_prj[4]), 
           expand = T) +
  scale_size(range = c(0.01, 2)) +
  # theme_void() +
  theme(
    plot.background = element_rect(fill = "white"),
    panel.background = element_rect(fill = "white"),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    panel.border = element_rect(colour = "black", fill=NA, size=.5),
    legend.key=element_blank()
  )

map



if(level == "site"){
} else if(level == "district"){ ## district (county)
  # ggsave(paste0("figures/", season, "_districts_tag_wash_ice_directX.png") plot=map, width=5, height = 8)
  ggsave(paste0("figures/", season, "_test_comb_euring_cring_directX.png"), plot=map, width=5.5, height = 7)
} else if(level == "country"){ ## Country
  ggsave(paste0("figures/", season, "_countries_tag_wash_ice_directX.png"), plot=map, width=5.5, height = 8)
  
}