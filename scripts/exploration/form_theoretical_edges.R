## forming edges based on model of movement probability between sites

pacman::p_load(dplyr, magrittr, sf, circular)


## save
allsites <- readRDS(
  "data/analysis/alldatatypes_allsites_poly_all.rds"
)

##
allsites$loc_num <- as.numeric(as.factor(allsites$SitRecID))  # (absolute) numeric

netdat %<>% left_join(allsites)

## only sites used during specific season
seasonsites <- allsites[allsites$SitRecID %in% unique(netdat$SitRecID), ] # netdat from network script (cleaned, season specific data)

## site centroids
site_cent <- st_centroid(seasonsites)

# site_cent <- site_cent[1:10, ]

#------------------------------------------------------------------------------
## Distance -- negative exponential function ----------------------------------
#------------------------------------------------------------------------------

## distance matrix among all sites --------------------------------------------
# distmtrx <- units::drop_units(st_distance(site_cent) / 1000) # km - spherical
distmtrx <- units::drop_units( # ellipsoidal dists
  lwgeom::st_geod_distance(site_cent, site_cent) / 1000) # km

rownames(distmtrx) <- site_cent$loc_num
colnames(distmtrx) <- site_cent$loc_num

## only upper 'triangle
upptri <- distmtrx[upper.tri(distmtrx)]


## ----------------------------------------------------------------------------
## probability of movement between nodes is function of distance
## negative exponential probability function, set using median displace. dist.
# prob <- exp(-(1/med_displ)*upptri)
# plot(upptri, prob)

## med_displ comes from 'step_lengths_btwn_site.R'
dist_prob <- exp(-(1/med_displ)*distmtrx)

## links above the maximum observed movement dist. have probability of 
# max_displ comes from 'step_lengths_btwn_site.R'
dist_prob[which(distmtrx > max_displ)] <- 0

## turn dist matrix into data.frame
allcombs <- expand.grid(rownames(distmtrx), colnames(distmtrx), stringsAsFactors = F)

allcombs$dist <- as.vector(distmtrx)       # distance
allcombs$dist_prob <- as.vector(dist_prob) # distance-based probability

allcombs %<>% 
  dplyr::rename(from = Var1, to = Var2)

## filter out self-links and zero prob. links
edgelist <- allcombs[-which(allcombs$from == allcombs$to), ]
edgelist %<>% filter(dist_prob != 0)


#------------------------------------------------------------------------------
## Direction -- von Mises distribution ----------------------------------------
#------------------------------------------------------------------------------

## remove links headed in wrong seasonal direction ----------------------------
## i.e. northward movements in fall and southwards in spring ------------------

# lwgeom::st_geod_azimuth() # sf friendly, but works on point sequences

## start-end points of all edges
startpnts <- left_join(
  edgelist, 
  mutate(site_cent, from = as.factor(loc_num)), by="from") %>% st_as_sf() %>% 
  as_Spatial()

endpnts <- left_join(
  edgelist, 
  mutate(site_cent, to = as.factor(loc_num)), by="to") %>% st_as_sf() %>% 
  as_Spatial()

## get azimuths (degrees)
edgelist$angle <- geosphere::bearing(startpnts, endpnts)

## rmv northward movements in fall and southwards in spring ------------------
# edgelist <- slice(edgelist, 1:10)

if(season == "spring"){
  northw <- which(edgelist$angle < 90 & edgelist$angle > -90)
  edgelist <- edgelist[northw, ]
} else if (season == "fall"){
  southw <- which(!(edgelist$angle < 90 & edgelist$angle > -90))
  edgelist <- edgelist[southw, ]
}


## Derive a probability of a movement also from it's direction ----------------

## shift negative degrees to positive to close 'wrapping' gap in dist.
angle_shift <- ifelse(edgelist$angle < 0, edgelist$angle + 360, edgelist$angle)

## empirical sample statistics
# angle_circ <- circular(edgelist$angle, units="degrees")
angle_circ <- circular(angle_shift, units="degrees")
circ_mean  <- mean.circular(angle_circ) # mu of von Mises
circ_sd <- sd.circular(angle_circ) # related to kappa of von Mises
circ_var <- var.circular(angle_circ)

## sample from von mises distribution to visualize theor. dist.
# sampled_vals <- rvonmises(
#   n=1000, 
#   mu=circ_mean, 
#   kappa=mle.vonmises(angle_circ)$kappa
# )
# hist(as.numeric(sampled_vals),20)

## probability density function - get probabilities for each edge-angle
dir_prob <- pvonmises(
  angle_circ, 
  mu=circ_mean, 
  kappa=mle.vonmises(angle_circ)$kappa
)
# plot(as.numeric(angle_circ), dir_prob) # slow

# plot(-180:180, 
#      pvonmises(
#        circular(-180:180, units="degrees"), 
#        mu=circ_mean, 
#        kappa=mle.vonmises(angle_circ)$kappa)
# )

## add angle probability by distance probability -------------------------

edgelist$dir_prob <- dir_prob
## average together the probabilities 
## NOT GREAT SINCE ONE IS EXPONENTIAL THERE ARE MULTIPLE ORDERS OF MAGNITUDE 
# DIFFERENCES IN PROBABILITIES

edgelist$prob <- (edgelist$dist_prob + edgelist$dir_prob)/2 ## average probs.

plot(edgelist$dist, edgelist$prob) # slow

edgelist$weight <- 1 - edgelist$prob # inverse (igraph: weight is a cost)

## Network --------------------------------------------------------------------

## undirected - reciprocal edges
netsf <- sfnetwork(nodelist, edgelist, node_key = "loc_num", directed = T, 
                   edges_as_lines = TRUE)
## directed - order matters
# netsf <- sfnetwork(nodelist, edgelist, node_key = "dist", directed = T, 
#                    edges_as_lines = TRUE)

## interactive map it
nodesf <- netsf %>% activate("nodes") %>% sf::st_as_sf()
edgesf <- netsf %>% activate("edges") %>% sf::st_as_sf()

# edgesf %>% mapview(zcol="angle")
# edgesf %>% slice(1800:2000) %>% mapview(zcol="angle")
edgesf %>% slice(2500:5000) %>% mapview(zcol="prob")
edgesf %>% slice(25000:50000) %>% mapview(zcol="dist")

mapview::mapview(edgesf, zcol="angle") + 
  mapview::mapview(nodesf) 
