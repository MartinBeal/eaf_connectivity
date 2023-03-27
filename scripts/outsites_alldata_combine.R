### ---------------------------------------------------------------------------
### Identify all outsites (hexcells) across all three datatypes ### -----------
## Combine hexcells outside IBAs into 'outsites' 
# Forms outsite layer to use in creating datatype-specific networks 
### ---------------------------------------------------------------------------

pacman::p_load(dplyr, stringr, ggplot2, sf, mapview, magrittr, lubridate, dggridR)

## fxn for splitting string into columns
source("C:/Users/Martim Bill/Documents/R/source_scripts/str2col.R")


# alldat <- readRDS("data/analysis/alldatatypes_100km_ibas.rds")
alldat <- readRDS("data/analysis/alldatatypes_100km_ibas_stpovrs.rds")

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

## Get the number of records in each equally-sized cell
binned_re <- outdat %>% group_by(cell) %>% summarise(count=n())

## Get the grid cell boundaries for cells which had obs
grid_re           <- dgcellstogrid(dggs, binned_re$cell)
colnames(grid_re)[1] <- "cell"
ccenters_re       <- dgSEQNUM_to_GEO(dggs, grid_re$cell)

# cbind.data.frame(lat=ccenters_re$lat_deg, lon=ccenters_re$lon_deg) %>%
#   sf::st_as_sf(
#     coords = c("lon", "lat"),
#     crs = 4326, agr = "constant") %>% mapview()

# mapview(grid_re, zcol = "count")

## add cell center coordinates
grid_df <- grid_re %>% st_drop_geometry() %>%
  bind_cols(latitude=ccenters_re$lat_deg, longitude=ccenters_re$lon_deg)

## Get the country that each cell center falls w/in ---------------------------

## country-eez union
# cntryeez <- st_read("C:/Users/Martim Bill/Documents/political_connectivity/spatial_data/shapefiles_EEZ_countries/union_countries_EEZs/EEZ_Land_v3_202030.shp")

## overlay points on polygons -------------------------------------------------
# grid_sf <- st_as_sf(
#   grid_df,
#   coords = c("longitude", "latitude"), 
#   crs = 4326, agr = "constant"
# )
# 
# ov <- sapply(
#   st_intersects(grid_sf, cntryeez),
#   function(x){
#     if(length(x) == 0){x <- 'none'}
#     return(x[1])
#   })
# 
# cntryeez$rowid <- 1:nrow(cntryeez)
# 
# ## combine site info with overlap result
# pntscntry <- left_join(data.frame(rowid = ov), st_drop_geometry(cntryeez))
# 
# ## combine overlap result back into bird locations w/ site info
# grid_sf <- bind_cols(grid_sf, pntscntry[,c("SOVEREIGN1", "ISO_SOV1")])
# 
# grid_sf <- rename(grid_sf, country = SOVEREIGN1, iso_a3 = ISO_SOV1)
# 
# ## Update the grid cells' properties to include country info
# grid_re %<>% right_join(st_drop_geometry(grid_sf))

## ----------------------------------------------------------------------------
## Group neighbouring hexcells into 'sites' -----------------------------------
## ----------------------------------------------------------------------------

## solution found here:
## https://gis.stackexchange.com/questions/323038/dissolve-only-overlapping-polygons-in-r-using-sf

## union intersecting polygons into combined polygon
parts <- st_cast(st_union(grid_re),"POLYGON")
plot(parts)

parts <- st_make_valid(parts)

## identify which 'part' each cell polygon intersects
cellgrp <- unlist(st_intersects(grid_re, parts))

## union cell polygons by their cellgroup membership + paste cell numbers
grid_re2 <- cbind(grid_re, cellgrp) %>%
  # group_by(cellgrp, country) %>%   ## split outsite by borders
  group_by(cellgrp) %>%            ## ignore country borders
  summarize(
    cell = paste(cell, collapse = ", ")
  )

mapview(grid_re2)

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

## give an identifier to share w/ IBA data ------------------------------------
grid_re2 %<>% mutate(
  cellgrp = paste0("cellgrp_", cellgrp),
  SitRecID  = cellgrp ## make the cell group the main site ID
) 

## Save layer of outsite (cell-group) polygons 

saveRDS(grid_re2, 
        "data/analysis/site_nodes/outsite_polygons_alldatatypes_all.rds")

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
  "data/analysis/combined/alldatatypes_ibas_stpovrs_outsites_all.rds")


# create reference table of sites w/ obs --------------------------------------
site_cent <- bind_rows(site_cent, out_cent)

site_summ <- alldat2 %>% 
  group_by(SitRecID) %>% summarise(site_poly = first(site_poly)) %>% 
  left_join(site_cent, by = c("SitRecID")) %>% st_as_sf()


## SAVE 
saveRDS(
  site_summ, 
  "data/analysis/site_nodes/alldatatypes_allsites_cent_all.rds")


## ----------------------------------------------------------------------------
## Combine outsite polygons w/ godwit IBA polygons ----------------------------
## polygon dataset (IBAs) 
iba <- st_read("data/geodata/ibas/EAF_btgo_IBA/eaf_btgo_iba.shp")

iba <- iba %<>% st_make_valid() # deal with polygon invalidity

useiba <- filter(iba, SitRecID %in% unique(site_summ$SitRecID))

useiba %<>% dplyr::select(SitRecID, Country, ISO3) %>% 
  mutate(SitRecID = as.character(SitRecID)) %>% 
  rename(country=Country)

# NEED TO PROJECT FIRST
outsites <- grid_re2
outsites$GISArea <- units::set_units(st_area(outsites), "km2")

## outsites 

outsites %<>% dplyr::select(SitRecID, GISArea, geometry)

## overlay site polygons on country polygons -----------------------------
cntryeez <- st_read("C:/Users/Martim Bill/Documents/political_connectivity/spatial_data/shapefiles_EEZ_countries/union_countries_EEZs/EEZ_Land_v3_202030.shp")

cntryeez$rowid <- as.character(1:nrow(cntryeez))

ov <- lapply(
  st_intersects(outsites, cntryeez),
  function(x){
    if(length(x) == 0){x <- 'none'}
    cntrys <- cntryeez[x, ]$ISO_SOV1
    return(paste(cntrys, collapse = "_"))
  })

## combine site info with overlap result
cntryiso <- do.call(rbind, ov)
cntrynames <- do.call(
  rbind,
  lapply(str_split(cntryiso, pattern = fixed("_")), function(x){
    countries <- unique(filter(cntryeez, ISO_SOV1 %in% x)$SOVEREIGN1)
    return(paste(countries, collapse = "_"))
  })
)

cntrys <- data.frame(country = cntrynames, ISO3 = cntryiso)

## combine overlap result back into bird locations w/ site info
outsites <- bind_cols(outsites, cntrys)

## combine w/ IBAs
allsites <- bind_rows(outsites, useiba)

## Save
saveRDS(allsites, 
        "data/analysis/site_nodes/allsite_polygons_alldatatypes_all.rds")
