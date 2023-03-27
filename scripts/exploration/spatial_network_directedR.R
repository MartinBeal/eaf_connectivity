### Create a spatial, DIRECTED network 

## 1: identify edges w/in season re-obs (tracking/color), weight by prop.
# of marked individuals
## 2: Form same edges for disconnected nodes as to it's closest connected neighbor
## 3: weight 'disconnected node' edges by ?empirical proportion and direction, 
# distance probability functions?

pacman::p_load(dplyr, igraph, stringr, tictoc, tidygraph, sfnetworks, ggplot2, 
               sf, mapview, magrittr, lubridate, netrankr)

## Load custom fxns
# fxn for splitting string into columns
source("C:/Users/Martim Bill/Documents/R/source_scripts/str2col.R")
# fxn for calculating graph-level centralization based on node centralization
source("C:/Users/Martim Bill/Documents/R/source_scripts/sfnet_globalmetrics.R")

## choose data subset to run --------------------------------------------------

## which network to create
# season <- "all"
# season <- "spring"
season <- "fall"

## Load site summary both IBAs and outsite centroids) ------------------------
site_summ <- readRDS(
  paste0("data/analysis/site_nodes/alldatatypes_allsites_cent_all.rds")
)

## Load combined data, with both IBA overlay info and hexgrid outsites -------
### SAVE 
alldat <- readRDS(
  paste0("data/analysis/combined/alldatatypes_ibas_stpovrs_outsites_all.rds"))

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
# netdat$loc_num <- as.numeric(as.factor(netdat$SitRecID))  # (absolute) numeric
# site_summ$loc_num <- as.numeric(as.factor(site_summ$SitRecID))  # (absolute) numeric
site_summ$loc_num <- site_summ$SitRecID # (absolute) numeric

netdat %<>% left_join(site_summ)


## Vertex list (all datatypes) ------------------------------------------------
n_id_total <- n_distinct(netdat$bird_id)

nodelist <- netdat %>% group_by(loc_num, SitRecID) %>% 
  summarise(
    n_id    = n_distinct(bird_id), # n IDs observed at site
    prop_id = n_id / n_id_total,   # prop. marked IDs
    n_obs   = n()
  ) #%>% right_join(n_obstype)

## Weight - proportion of marked population using node
# nodelist$weight <- nodelist$prop_id

nodelist <- nodelist %>% 
  left_join(site_summ) %>% 
  sf::st_as_sf()


## Edge list ------------------------------------------------------------------
## Defining edges based on w/in month obs-reobs -------------------------------

oneid_list <- split(netdat, netdat$bird_id)

oneid_list <- lapply(
  seq_along(oneid_list), 
  function(x){
    print(x)
    one <- oneid_list[[x]]
    
    from <- one[1:nrow(one)-1,]
    to   <- one[2:nrow(one),]
    
    xx <- data.frame(
      from  = from$loc_num,
      to    = to$loc_num,
      from_sid  = from$SitRecID,
      to_sid   = to$SitRecID,
      bird_id  = one$bird_id[1],
      datatype = one$datatype[1],
      tdiff = difftime(to$timestamp, from$timestamp, units = "day")
    )
    
    return(xx)
    
  }
)

## each obs is a new link
obsfull <- data.table::rbindlist(oneid_list)

## time span between re-obs (i.e., forming the link)
ggplot() +
  geom_histogram(data=obsfull, aes(x=as.numeric(tdiff))) +
  facet_wrap(~datatype, scales = "free_x")

## retain only steps of < 32 days time difference 
max_time <- 32
obsred <- filter(obsfull, tdiff < max_time)

## remove self connections
obsred <- obsred[-which(obsred$from == obsred$to), ]

## short-term time span between re-obs (i.e., forming the link)
ggplot() + 
  geom_histogram(data=obsred, aes(x=as.numeric(tdiff))) + 
  facet_wrap(~datatype, scales = "free_x") + theme_bw() +
  xlab("Days between observations")

## SAVE
# ggsave(
#   paste0("figures/step_under_32days_datatype_", season, ".png"), 
#   width = 9, height = 3)


## Population-level edges only (i.e. no overlapping links) --------------------
# identify identical connections ignoring order (undirected)
obsred$sortcomb <- sapply(
  seq_along(obsred$from), function(f){
    rw <- obsred[f,]
    # sortcomb <- paste(sort(c(rw$from, rw$to)), collapse = " ")
    sortcomb <- paste(rw$from, rw$to)
    return(sortcomb)
  })

### 'flat' pop-level edges - only one for site-site connection (directional)
if(any(duplicated(obsred$sortcomb))){
  
  ## summarize by datatype
  dtype_summ <- obsred %>% group_by(from, to, sortcomb, datatype) %>% 
    summarise(
      n_bird = n_distinct(bird_id)
    )
  
  popflat <- dtype_summ %>% group_by(from, to, sortcomb) %>% 
    summarise(
      n_birds   = sum(n_bird),
      datatypes = paste(unique(datatype), collapse = " "),
    )
  
  ## only tracking (b/c n tracked birds forming link least biased measure)
  # trx_summ <- filter(dtype_summ, datatype == "trax")
  #   ## add n
  # popflat %<>% 
  #   left_join(trx_summ[, c("sortcomb", "n_bird")]) %>% 
  #   rename(n_bird_trax = n_bird)
}

edgelist <- popflat %>% dplyr::rename(link_id = sortcomb) %>% 
  as_tibble()

## make location codes character to allow for 'missing' links 
nodelist$loc_num <- as.character(nodelist$loc_num)

edgelist %<>% mutate(
  from = as.character(from),
  to   = as.character(to),
  from_sid = from, # maintaining global site codes
  to_sid = to      # maintaining global site codes
)

## Sum in-birds and out-birds for each node from edge data --------------------

from_birds <- edgelist %>% group_by(from_sid) %>% 
  summarise(
    n_id_out = sum(n_birds)
  )

in_birds <- edgelist %>% group_by(to_sid) %>% 
  summarise(
    n_id_in = sum(n_birds)
  )

inouts <- left_join(from_birds, in_birds, by = c("from_sid" = "to_sid")) %>% 
  rename(SitRecID = from_sid)

## add to nodelist
nodelist %<>% left_join(inouts)

## add to edgelist
edgelist %<>% left_join(from_birds) %>% 
  rename(
    n_id_out_from = n_id_out
  ) %>% 
  mutate(
    prop_id_out_from = n_birds / n_id_out_from
  )

## Create sfnetwork -----------------------------------------------------------

netsf <- sfnetwork(nodelist, edgelist, directed = TRUE, 
                   edges_as_lines = TRUE)

## ego metrics (node/edge-level)
netsf %<>%
  activate(nodes) %>%
  mutate(
    degree_tot      = centrality_degree(loops = FALSE), # no self-connections
    degree_in       = centrality_degree(loops = FALSE, mode = "in"), # no self-connections
    degree_out      = centrality_degree(loops = FALSE, mode = "out"), # no self-connections
    # degree_norm = degree / n_distinct(loc_num), # normalized 0-1 (prop of sites)
    degree_tot_rank = dense_rank(desc(degree_tot)),
    degree_in_rank  = dense_rank(desc(degree_in)),
    degree_out_rank = dense_rank(desc(degree_out)),
    # strength    = centrality_degree(loops = FALSE, weights = weight), # sum of edge weights
    between     = centrality_betweenness(directed = TRUE, weights = NULL), # node betweenness
    # between_w   = centrality_betweenness(directed = FALSE, weights = weight), # node betweenness
    # between_norm = centrality_betweenness(directed = TRUE, normalized = T), # normalized
    btwn_rank   = dense_rank(desc(between)),
    # btwn_w_rank   = dense_rank(desc(between_w))
  )


nodesf <- netsf %>% activate("nodes") %>% sf::st_as_sf()
edgesf <- netsf %>% activate("edges") %>% sf::st_as_sf()

# mapview::mapview(edgesf, zcol="n_bird_trax") + 
#   mapview::mapview(nodesf) 

# mapview::mapview(edgesf) + 
#   mapview::mapview(nodesf, zcol = "btwn_rank") 
# 
# mapview::mapview(nodesf, zcol = "n_id_out") 
# 
# mapview::mapview(edgesf, zcol = "prop_id_out_from") 

## disconnected nodes
disc_nod_ids <- unique(nodesf$SitRecID)[which(
  !unique(nodesf$SitRecID) %in% as.character(unique(c(edgesf$from_sid, edgesf$to_sid)))
  )]

disc_nod <- filter(nodesf, loc_num %in% disc_nod_ids)

disc_nod %>% mapview() # all
filter(disc_nod, n_id > 1) %>% mapview() # sites w/ more than 1 tracked ID
filter(disc_nod, n_obs > 1) %>% mapview() # sites w/ more than 1 tracking period

## SAVE ##

saveRDS(
  netsf, 
  paste0("data/analysis/networks/directed/", max_time,"d_", season, "_iba_hex_10km.rds")
)
