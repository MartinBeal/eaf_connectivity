### overlay bird locations on vector dataset of 'sites' (IBAs, Ramsar, PAs)
## and find closest site w/in a threshold distance
# Also quantify # consec. days bird is at same site
## ----------------------------------------------------------------------------

pacman::p_load(dplyr, igraph, stringr, tictoc, ggplot2, data.table,
               sf, mapview, magrittr, lubridate)


## tracking data
# alldat <- readRDS("data/analysis/tracking/PTT_GPS_mconn_12h.rds")
alldat <- readRDS("data/analysis/tracking/PTT_GPS_mconn_12h_no0.rds")

alldat %<>% rename(bird_id = id)

## polygon dataset (IBAs)
# iba <- st_read("data/geodata/ibas/EAF_btgo_IBA/eaf_btgo_iba.shp")
iba <- st_read("data/geodata/ibas/Africa_Europe_IBA/Africa_Europe_IBA.shp")

iba %<>% st_make_valid() # deal with polygon invalidity

## get centroids of sites for later use as nodes in the network
iba_cent <- iba %>% st_centroid()

saveRDS(iba_cent, "data/geodata/ibas/Africa_Europe_IBA/Africa_Europe_IBA_centroids.shp")

## small subset 
# alldat <- head(alldat, 1000)

alldat_sf <- alldat %>% st_as_sf(
  coords = c("longitude", "latitude"), 
  crs = 4326)

## overlay points on polygons -------------------------------------------------
ov <- sapply(
  st_intersects(alldat_sf, iba),
  function(x){
    if(length(x) == 0){x <- 'none'}
    return(x[1])
  })

iba$rowid <- as.character(1:nrow(iba))

## combine site info with overlap result
pntsibas <- left_join(data.frame(rowid = ov), st_drop_geometry(iba))

## combine overlap result back into bird locations w/ site info
alldat <- bind_cols(alldat, pntsibas[,c("SitRecID", "IntName")])

## name points outside polygons as 'none'
alldat %<>% mutate(
  SitRecID = ifelse(is.na(SitRecID), "none", SitRecID),
  IntName = ifelse(is.na(IntName), "none", IntName)
)


# alldat_sf <- alldat %>% st_as_sf(
#   coords = c("longitude", "latitude"), 
#   crs = 4326)
# mapview(subset(iba, SitRecID %in% unique(alldat$SitRecID))) + mapview(alldat_sf, zcol="IntName")


## get nearest polygon(s) for points outside ----------------------------------

out <- subset(alldat, SitRecID == "none")
# out <- head(out, 1000)

out_sf <- out %>% st_as_sf(
  coords = c("longitude", "latitude"), 
  crs = 4326)

tictoc::tic()
ibainrange <- st_is_within_distance(out_sf, iba, 10000) ## m from near. poly
tictoc::toc()

## classify result (none, iba row id, multiple)
whichpoly <- sapply(
  ibainrange,
  function(x){
    if(length(x) == 0)
    {x <- 'none'} else if (length(x) > 1)
    {x <- 'multiple'
    }
    return(x[1])
  })

## get rowid of ibas for points falling in just one iba (keep 'nones')
singles_rid <- as.character(whichpoly[which(whichpoly != "multiple")])
singles <- out[which(whichpoly != "multiple"), ] # filter from bird data

## get iba metadata for ibas that points fall within
single_ibas <- left_join(
  data.frame(rowid = singles_rid), 
  st_drop_geometry(iba)[, c("rowid", "SitRecID", "IntName")]) %>% 
  mutate(
    SitRecID = as.character(SitRecID))

## combine overlap result back into bird locations w/ site info (#2)
singles <- bind_cols(
  subset(singles, select=-c(SitRecID, IntName)),
  single_ibas[,c("SitRecID", "IntName")])


## points overlapping multiple polygons ---------------------------------------
multis   <- out[which(whichpoly == "multiple"), ]
multi_sf <- out_sf[which(whichpoly == "multiple"), ]

## get nearest polygon for those obs w/ multiple within search distance
tictoc::tic()
neariba <- st_nearest_feature(multi_sf, iba)
tictoc::toc()

## combine site info with overlap result
multis_ibas <- left_join(
  data.frame(rowid = as.character(neariba)), 
  st_drop_geometry(iba)[, c("rowid", "SitRecID", "IntName")]) %>% 
  mutate(
    SitRecID = as.character(SitRecID))

## combine overlap result back into bird locations w/ site info (#3)
multis <- bind_cols(
  subset(multis, select=-c(SitRecID, IntName)),
  multis_ibas[,c("SitRecID", "IntName")])


## recombine all data and re-arrange ------------------------------------------
## remove outs from all points
all_noout <- subset(alldat, SitRecID != "none")

alldat2 <- bind_rows(all_noout, singles, multis) %>% 
  dplyr::arrange(bird_id, timestamp)

## again name points outside polygons as 'none'
alldat2 %<>% mutate(
  SitRecID = ifelse(is.na(SitRecID), "none", SitRecID),
  IntName = ifelse(is.na(IntName), "none", IntName)
)

### do some sense-checks
## check point falling in/around a site
# mapview(subset(iba, IntName == "Tejo estuary")) + (alldat2 %>% filter(IntName == "Tejo estuary") %>% st_as_sf(
#   coords = c("longitude", "latitude"),
#   crs = 4326) %>% mapview())

# ## points not considered in/around a site 
# mapview(subset(iba, IntName == "Tejo estuary")) + (alldat2 %>% filter(IntName == "none") %>% st_as_sf(
#   coords = c("longitude", "latitude"),
#   crs = 4326) %>% mapview())


## Calculate number of consecutive days spent at a site -----------------------
alldat2 <- as.data.table(alldat2)

## id consecutive obs at a site, per bird (data.table way)
alldat2[, samesite := data.table::rleid(SitRecID),
       by = bird_id]
alldat2$samesite <- ifelse(is.na(alldat2$SitRecID), NA, alldat2$samesite) # NA where site NA

## number of consecutive days at a site
alldat2 <- alldat2 %>%
  group_by(bird_id, samesite) %>%
  summarise(n_day = n_distinct(yday(timestamp))) %>% 
  left_join(alldat2)

## remove points outside polygons and points w/in from visits <= 1 day
# alldat2 %>% filter(n_day > 1 & SitRecID != "none") %>% st_as_sf(
#   coords = c("longitude", "latitude"),
#   crs = 4326) %>% mapview()

# alldat2 %>% filter(bird_id == "R\xfcschke") %>%
#   st_as_sf(coords = c("longitude", "latitude"),
#            crs = 4326, agr = "constant") %>% mapview()

## rename
alldat2 %<>% rename(site_poly = IntName) ## IBAs only

## SAVE -----------------------------------------------------------------------

# saveRDS(alldat2, "data/analysis/tracking/PTT_GPS_mconn_12h_ibas.rds")
saveRDS(alldat2, "data/analysis/tracking/PTT_GPS_mconn_12h_no0_ibas.rds")


## get distances 
# tictoc::tic()
# dist <- st_distance(out_sf, iba[neariba,], by_element=TRUE)
# tictoc::toc()