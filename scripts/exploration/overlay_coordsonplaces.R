## Convert all coordinates to the lowest accuracy in the data (e.g. county)
## Overlay coordinates on map of sites/areas/counties/countries and extract

pacman::p_load(sf)

# alldat <- readRDS("data/analysis/cr_tag_ice_wash_merge_newcoords.rds")  # c-ringing
alldat <- readRDS("data/analysis/test_comb_euring_cring_newcoords.rds") # combined euring/cring

## vectorized fxn to count number of decimal place (https://stackoverflow.com/questions/5173692/how-to-return-number-of-decimal-places-in-r)
decimalplaces <- function(x) {
  ifelse(abs(x - round(x)) > .Machine$double.eps^0.5,
         nchar(sub('^\\d+\\.', '', sub('0+$', '', as.character(x)))),
         0)
}

dplace <- decimalplaces(alldat$longitude)
dplace <- decimalplaces(alldat$latitude)

alldat$latitude_r  <- round(alldat$latitude, 3)
alldat$longitude_r <- round(alldat$longitude, 3)

##

places <- alldat %>% group_by(x, y) %>% # x, y is address
  summarise(
    longitude = first(longitude_r),
    latitude = first(latitude_r)
  )

# write.csv(places, "data/analysis/counties.csv", row.names = F)

## map locations
coo_sf <- st_as_sf(
  alldat, coords = c("longitude", "latitude"),
  crs = 4326, agr = "constant")
# mapview::mapview(coo_sf)


### overlay obs coords on districts layer -------------------------------------

distr <- st_read("data/geodata/EAF_districts_21JUL22/EAF_districts_21JUL22.shp")

distr <- st_make_valid(distr)

# ovs <- st_intersects(coo_sf, distr, sparse = T)
# ovs <- st_intersects(coo_sf, distr, sparse = T)
tictoc::tic()
ovs <- st_nearest_feature(coo_sf, distr)
tictoc::toc()

dists <- do.call(
  rbind,
  lapply(ovs, function(x){
    if(length(x) == 0){
      xx <- NA
      } else {
        xx <- x
      }
    if(!is.na(xx)){
      onedist <- distr[xx, ]
      if(is.na(onedist$NAME_2)){ name <- onedist$NAME_1 } else{
        name <- onedist$NAME_2
      } 
      return(name)
    } else {  return(NA)}
    })
)

alldat$dist <- as.vector(dists)

## SAVE ## 
# saveRDS(alldat, "data/analysis/cr_tag_ice_wash_merge_districts.rds")
saveRDS(alldat, "data/analysis/test_comb_euring_cring_districts.rds")
