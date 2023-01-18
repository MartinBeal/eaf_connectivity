## Group neighbouring hexcells into 'sites' 

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
  bind_cols(latitude=ccenters_re$lat_deg, longitude=ccenters_re$lon_deg)


### Overlay relocations on cell group polygons --------------------------------
outdat %<>% sf::st_as_sf(
  coords = c("longitude", "latitude"),
  crs = 4326, agr = "constant")

ov <- sapply(
  st_intersects(outdat, grid_re2),
  function(x){
    if(length(x) == 0){x <- 'none'}
    return(x[1])
  })

grid_re2$rowid <- 1:nrow(grid_re2)

## combine site info with overlap result
pntscellgrps <- left_join(data.frame(rowid = ov), st_drop_geometry(grid_re2))

## combine overlap result back into bird locations w/ cellgroup info
outdat <- bind_cols(outdat, pntscellgrps[,c("cellgrp", "cell")]) %>% 
  mutate(
    cellgrp = paste0("cellgrp_", cellgrp),
    SitRecID = cellgrp ## make the cell group the main site ID
  )

## create reference table of sites w/ obs --------------------------------------
## Out site centroids:
out_cent <- ccenters_re %>% 
  mutate(
    cellgrp = paste0("cellgrp_", cellgrp),
    SitRecID = cellgrp ## make the cell group the main site ID
  ) %>% group_by(SitRecID, cellgrp) %>% 
  st_as_sf(coords = c("longitude", "latitude"), 
           crs = 4326, agr = "constant") %>% 
  summarise() %>% 
  dplyr::select(SitRecID, geometry)


## centroids of site polygons
site_cent <- readRDS( 
  "data/geodata/ibas/Africa_Europe_IBA/Africa_Europe_IBA_centroids.shp") 
site_cent %<>% 
  mutate(SitRecID = as.character(SitRecID)) %>% 
  dplyr::select(SitRecID, geometry)

## re-merge in and out data
netdat <- outdat %>% st_drop_geometry() %>% bind_rows(indat) %>%
  arrange(bird_id, timestamp)

# create reference table of sites w/ obs --------------------------------------
site_cent <- bind_rows(site_cent, out_cent)

site_summ <- netdat %>% 
  group_by(loc_num, site_poly, SitRecID) %>% summarise() %>% 
  left_join(site_cent, by = c("SitRecID")) %>% st_as_sf()
