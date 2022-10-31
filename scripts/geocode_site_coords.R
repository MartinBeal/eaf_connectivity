### Get coordinates for addresses (sites, counties, countries etc.)

pacman::p_load(tidygeocoder, stringr, dplyr)

# alldat <- readRDS("data/analysis/cr_tag_ice_wash_merge.rds")  # cringing
alldat <- readRDS("data/analysis/test_comb_euring_cring.rds") # combined euring/cring

# fxn for splitting string into columns
source("C:/Users/Martim Bill/Documents/R/source_scripts/str2col.R")

alldat$rownum <- seq_len(nrow(alldat))

## Use geocoding to get coordinates for obs missing coords ---------------------
origrws <- which(with(alldat, is.na(latitude) | is.na(longitude)))
nocoord <- alldat[origrws, ]

## set up which two columns will be in address 
nocoord$x <- with(nocoord, ifelse(is.na(county), area, county))
## use county where country missing
rws <- which(is.na(nocoord$country))
nocoord$y <- nocoord$country
## use 'area' where county missing
nocoord[rws,]$x <- nocoord[rws,]$area
nocoord[rws,]$y <- nocoord[rws,]$county
## bespoke fixes
nocoord$x <- ifelse(nocoord$x == "Edinburghshire", "Midlothian", nocoord$x)

## remove any leftovers w/ no info to use
nocoord <- nocoord[!is.na(nocoord$x),]

## summary of unique addresses
counties <- nocoord %>% select(x, y) %>% 
  group_by(x) %>% summarise(y=first(y))

## geocode coordinates from addresses -----------------------------------------
# councoord <- geo(address = counties$x)
# councoord2 <- geo(county = counties$x, country = counties$y)
# councoord3 <- geo(address = paste(counties$x, counties$y, sep=", "))
councoord4 <- geo(address = paste(counties$x, counties$y, sep=", "), method = "arcgis")

## map locations
# coo_sf <- st_as_sf(councoord4, coords = c("long", "lat"),
#          crs = 4326, agr = "constant")
# mapview::mapview(coo_sf)

## split address and clean up
newcoords <- cbind(
  str2col(
    councoord4$address, pattern=", ", 
    cols=1:2, colnames=c("x", "y")),
  councoord4
) %>% 
  dplyr::select(-address)

## add new coords back in
wcoords <- left_join(nocoord, newcoords) %>% 
  select(-latitude, -longitude) %>% 
  rename(latitude = lat, longitude = long) #%>% 
  # select(-x, -y)

## remove old rows and replace w/ same data including new coords
subdat <- alldat[-origrws,]
subdat %<>% mutate(x = county, y = country)
alldat_out <- rbind.data.frame(subdat, wcoords) %>% arrange(rownum)

## SAVE ## 
# saveRDS(alldat_out, "data/analysis/cr_tag_ice_wash_merge_newcoords.rds")
saveRDS(alldat_out, "data/analysis/test_comb_euring_cring_newcoords.rds")
