### ---------------------------------------------------------------------------
### Identify all outsites (hexcells) across all three datatypes ### -----------
## Combine hexcells outside IBAs into 'outsites' 
# Forms outsite layer to use in creating datatype-specific networks 
### ---------------------------------------------------------------------------

pacman::p_load(dplyr, stringr, ggplot2, sf, mapview, magrittr, lubridate, dggridR)

## fxn for splitting string into columns
source("C:/Users/Martim Bill/Documents/R/source_scripts/str2col.R")

## ring relocations overlaid on polygon layer
coldat <- readRDS("data/analysis/ringing/cring_merge_no7dayreobs_ibas.rds")
coldat %<>% rename(timestamp = date) %>% mutate(datatype="color")
## metal ring locations overlaid on polygon layer
metdat <- readRDS("data/analysis/ringing/euring_metal_ibas.rds")
metdat %<>% rename(timestamp = date) %>% mutate(datatype="metal")
## tracking locations overlaid on polygon layer
trxdat <- readRDS("data/analysis/tracking/PTT_GPS_mconn_12h_no0_stpovrs_ibas.rds")
trxdat %<>% mutate(datatype="trax")


## combine datatype sets ------------------------------------------------------
alldat <- bind_rows(
  coldat[, c("bird_id", "timestamp", "latitude", "longitude", "datatype", "site_poly", "SitRecID")],
  metdat[, c("bird_id", "timestamp", "latitude", "longitude", "datatype", "site_poly", "SitRecID")],
  trxdat[, c("bird_id", "timestamp", "latitude", "longitude", "datatype", "site_poly", "SitRecID")]
)

## Use only locations falling outside any existing site polygon
outdat <- alldat %>% filter(site_poly == "none")
indat  <- alldat %>% filter(site_poly != "none")

## ----------------------------------------------------------------------------
## Hexgrid --------------------------------------------------------------------
## ----------------------------------------------------------------------------

## Construct global grid with cells approximately 1000 miles across
# dggs <- dgconstruct(spacing=8, resround='down')
dggs <- dgconstruct(spacing=10, resround='down')

## create hexgrid from cr obs
outdat$cell <- dgGEO_to_SEQNUM(
  dggs, outdat$longitude, outdat$latitude)$seqnum

#Get the number of records in each equally-sized cell
binned_re <- outdat %>% group_by(cell) %>% summarise(count=n())

#Get the grid cell boundaries for cells which had obs
grid_re           <- dgcellstogrid(dggs, binned_re$cell)
colnames(grid_re)[1] <- "cell"
ccenters_re       <- dgSEQNUM_to_GEO(dggs, grid_re$cell)

# cbind.data.frame(lat=ccenters_re$lat_deg, lon=ccenters_re$lon_deg) %>%
#   sf::st_as_sf(
#     coords = c("lon", "lat"),
#     crs = 4326, agr = "constant") %>% mapview()

## Update the grid cells' properties to include the number of obs in each cell
grid_re <- merge(grid_re, binned_re, by.x="cell", by.y="cell")

# mapview(grid_re, zcol = "count")

## add cell center coordinates
grid_df <- grid_re %>% st_drop_geometry() %>%
  bind_cols(latitude=ccenters_re$lat_deg, longitude=ccenters_re$lon_deg)


## ----------------------------------------------------------------------------
## Group neighbouring hexcells into 'sites' -----------------------------------
## ----------------------------------------------------------------------------

## solution found here:
## https://gis.stackexchange.com/questions/323038/dissolve-only-overlapping-polygons-in-r-using-sf

## union intersecting polygons into combined polygon
parts <- st_cast(st_union(grid_re),"POLYGON")
plot(parts)

## identify which 'part' each cell polygon intersects
cellgrp <- unlist(st_intersects(grid_re, parts))

## union cell polygons by their cellgroup membership + paste cell numbers
grid_re2 <- cbind(grid_re, cellgrp) %>%
  group_by(cellgrp) %>%
  summarize(
    cell = paste(cell, collapse = ", ")
  )

# mapview(grid_re2)

## centroid of cell group polygon
ccenters_re <- st_centroid(grid_re2) %>% 
  dplyr::mutate(lon_deg = sf::st_coordinates(.)[,1],
                lat_deg = sf::st_coordinates(.)[,2])

## interactive map of cell centers
# cbind.data.frame(lat=ccenters_re$lat_deg, lon=ccenters_re$lon_deg) %>%
#   sf::st_as_sf(
#     coords = c("lon", "lat"),
#     crs = 4326, agr = "constant") %>% mapview()

### Update the grid cells' properties to include the number of obs in each cell
# grid_re2 <- merge(grid_re2, binned_re, by.x="cell", by.y="cell")
# mapview(grid_re, zcol = "count")
# mapview(grid_re, zcol = "n_birds")

## add cell center coordinates
grid_re2 <- grid_re2 %>% 
  bind_cols(
    latitude = ccenters_re$lat_deg, longitude = ccenters_re$lon_deg)


### Overlay relocations on cell group polygons --------------------------------
outdat_sf <- outdat %>% sf::st_as_sf(
  coords = c("longitude", "latitude"),
  crs = 4326, agr = "constant") %>% 
  bind_cols(outdat[,c("longitude", "latitude")])

ov <- sapply(
  st_intersects(outdat_sf, grid_re2),
  function(x){
    if(length(x) == 0){x <- 'none'}
    return(x[1])
  })

grid_re2$rowid <- 1:nrow(grid_re2)

## combine site info with overlap result
pntscellgrps <- left_join(data.frame(rowid = ov), st_drop_geometry(grid_re2))

## combine overlap result back into bird locations w/ cellgroup info
outdat_sf <- bind_cols(outdat_sf, pntscellgrps[, c("rowid", "cellgrp")]) %>% 
  mutate(
    cellgrp = paste0("cellgrp_", cellgrp),
    SitRecID = cellgrp ## make the cell group the main site ID
  )

## Save layer of outsite (cell-group) polygons 

saveRDS(grid_re2, 
        paste0("data/analysis/site_nodes/outsite_polygons_alldatatypes_", season, ".rds"))

## Recombine out and in data w/ site info for all data ------------------------

## Out site centroids:
out_cent <- ccenters_re %>% 
  mutate(
    cellgrp = paste0("cellgrp_", cellgrp),
    # site_poly = cellgrp, ## make the cell group the main site ID
    SitRecID  = cellgrp ## make the cell group the main site ID
  ) %>% group_by(SitRecID, cellgrp) %>% 
  st_as_sf(coords = c("longitude", "latitude"), 
           crs = 4326, agr = "constant") %>% 
  summarise() %>% 
  dplyr::select(SitRecID, geometry)

## centroids of IBA polygons
site_cent <- readRDS( 
  "data/geodata/ibas/Africa_Europe_IBA/Africa_Europe_IBA_centroids.shp") 

site_cent %<>% 
  mutate(SitRecID = as.character(SitRecID)) %>% 
  dplyr::select(SitRecID, geometry)#%>% 
  # rename(site_poly = IntName)

## re-merge in and out data
alldat2 <- outdat_sf %>% st_drop_geometry() %>% bind_rows(indat) %>%
  arrange(bird_id, timestamp)

### SAVE 
saveRDS(
  alldat2, 
  paste0("data/analysis/combined/alldatatypes_ibas_outsites_", season, ".rds"))


# create reference table of sites w/ obs --------------------------------------
site_cent <- bind_rows(site_cent, out_cent)

site_summ <- alldat2 %>% 
  group_by(SitRecID) %>% summarise(site_poly = first(site_poly)) %>% 
  left_join(site_cent, by = c("SitRecID")) %>% st_as_sf()


## SAVE 
saveRDS(
  site_summ, 
  paste0("data/analysis/site_nodes/alldatatypes_allsites_cent_", season, ".rds"))
