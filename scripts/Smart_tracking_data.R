
pacman::p_load(dplyr, sf, mapview, stringr, data.table, magrittr, amt)


## Dutch team tracking data 
raw  <- read.csv("data/tracking/team_smart/geo_trax.csv")
meta <- readxl::read_xlsx("data/tracking/team_smart/UK Project Godwit data.xlsx", sheet=1)

raw %<>% mutate(
  latitude  = as.numeric(lat),
  longitude = as.numeric(lon)
)

raw %<>% filter(!is.na(latitude))
raw %<>% filter(!is.na(longitude))

raw %>% sf::st_as_sf(
  coords = c("lon", "lat"), 
  crs = 4326, agr = "constant") %>%
  mapview(zcol="cr")
