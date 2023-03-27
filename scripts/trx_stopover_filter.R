### For tracking data, rmv non-stopover data
## bin day by hex cells and omit periods when bird didn't spend >= 2 d in a a site
# or cell

pacman::p_load(dplyr, stringr, ggplot2, sf, mapview, magrittr, lubridate, dggridR)

## fxn for splitting string into columns
source("C:/Users/Martim Bill/Documents/R/source_scripts/str2col.R")

### Season
season <- "all"
# season <- "spring"
# season <- "fall"

## datatypes combined 
alldat <- readRDS("data/analysis/alldatatypes_100km_ibas.rds")
trxdat <- filter(alldat, datatype == "trax")

## OLD
# trxdat <- readRDS("data/analysis/tracking/PTT_GPS_mconn_12h_no0_ibas.rds")

### Separate networks for fall and spring migration ---------------------------
## Spring: January 1 - June 30th, Fall: June 24th - January 31st

if(season == "spring"){
  trxdat <- subset(trxdat, month(trxdat$timestamp) %in% c(1:6))
} else if(season == "fall"){
  doy <- lubridate::yday(trxdat$timestamp) # June 23/24: 175
  trxdat <- subset(
    alldat, 
    month(trxdat$timestamp) %in% c(1, 7:12) | doy %in% c(175:181)
  )
}

## Use only locations falling outside any existing site polygon
nones <- trxdat %>% 
  filter(site_poly == "none")


## Hexgrid --------------------------------------------------------------------

## Construct global grid with cells approximately 1000 miles across
# dggs <- dgconstruct(spacing=8, resround='down')
dggs <- dgconstruct(spacing=10, resround='down')

## create hexgrid from cr obs
nones$cell <- dgGEO_to_SEQNUM(
  dggs, nones$longitude, nones$latitude)$seqnum

#Get the number of records in each equally-sized cell
binned_re <- nones %>% group_by(cell) %>% summarise(count=n())

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
grid_re <- grid_re %>% st_drop_geometry() %>%
  bind_cols(latitude=ccenters_re$lat_deg, longitude=ccenters_re$lon_deg)


## Calculate number of consecutive days spent at a site -----------------------
nones <- as.data.table(nones)
nones[, c("samesite","n_day"):=NULL]

## id consecutive obs at a site, per bird (data.table way)
nones[, samesite := data.table::rleid(cell),
       by = bird_id]
nones$samesite <- ifelse(is.na(nones$cell), NA, nones$samesite) # NA where site NA

nones %<>% 
  group_by(bird_id, samesite) %>%
  summarise(n_day = n_distinct(yday(timestamp))) %>%
  left_join(nones)

## number of consecutive days at a site
trxdat2 <- trxdat %>% 
  filter(site_poly != "none") %>% 
  bind_rows(nones) %>% 
  arrange(bird_id, timestamp)


## remove periods outside >= 48h stay at a site/cell
trxdat2 <- filter(trxdat2, n_day >= 2)


alldat2 <- bind_rows(trxdat2, filter(alldat, datatype != "trax")) 

## SAVE 

# saveRDS(alldat2, "data/analysis/tracking/PTT_GPS_mconn_12h_no0_stpovrs_ibas.rds")
saveRDS(alldat2, "data/analysis/alldatatypes_100km_ibas_stpovrs.rds")
