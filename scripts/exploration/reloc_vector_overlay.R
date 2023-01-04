## overlay bird locations on vector dataset of 'sites' (IBAs, Ramsar, PAs)
# and find closest site w/in a threshold distance

pacman::p_load(dplyr, igraph, stringr, tictoc, ggplot2, 
               sf, mapview, magrittr, lubridate)

datatype <- "metal"
# datatype <- "color"

if(datatype == "color"){
  ## ring relocations overlaid on polygon layer
  alldat <- readRDS("data/analysis/ringing/cring_merge_no7dayreobs.rds")
} else if (datatype == "metal"){
  ## tracking locations overlaid on polygon layer
  alldat <- readRDS("data/analysis/ringing/euring_filtbycr.rds") ## IBAs
}

## polygon dataset (IBAs)
iba <- st_read("data/geodata/ibas/EAF_btgo_IBA/eaf_btgo_iba.shp")

iba %<>% st_make_valid() # deal with polygon invalidity

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


## multiple polygons ---------------------------------------------------------
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
  arrange(bird_id, date)

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

## SAVE -----------------------------------------------------------------------

if(datatype == "color"){
  ## color resightings
  saveRDS(alldat2, "data/analysis/ringing/cring_merge_no7dayreobs_ibas.rds")
  } else if (datatype == "metal"){
  ## tracking locations overlaid on polygon layer
    saveRDS(alldat2, "data/analysis/ringing/euring_metal_ibas.rds")
}

# saveRDS(alldat2, "data/analysis/ringing/comb_euring_cring_no7dayreobs_ibas.rds")


## get distances 
# tictoc::tic()
# dist <- st_distance(out_sf, iba[neariba,], by_element=TRUE)
# tictoc::toc()