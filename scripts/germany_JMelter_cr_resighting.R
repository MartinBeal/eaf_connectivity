### Load German Schleswig-Holstein capture, recapture and resighting database - clean

## load packages
pacman::p_load(dplyr, magrittr, lubridate, stringr)

gdu <- readxl::read_xlsx('data/color/JMelter_germany_dummer/Sightings_Duemmer-birds-enroute.xlsx')
gdu <- readxl::read_xlsx('data/color/JMelter_germany_dummer/Data-Duemmer-new.xlsx')

# 
# ## translate column names 
# gdu %<>% dplyr::rename(
#   combo = `colour-ring /German-Code`, metal = `Ring-Nr`, 
#   latitude = `geografische Breite (Lat.)`, 
#   longitude = `geografische Länge (Long.)`, site = Site
# )

## translate column names 
gdu %<>% dplyr::rename(
  metal = `Ring-1`, 
  date = `date sighting`,
  latitude = `geografische Breite`, 
  longitude = `geografische Länge`, age = "age1",
) %>% mutate(
  date = dmy(date, tz="UTC"),
  latitude = as.numeric(str_replace(latitude, ",", ".")),
  longitude = as.numeric(str_replace(longitude, ",", "."))
  )

## create ringing observation for each individual based on data and Dummer lat/lon
ringrecs <- gdu %>% group_by(metal) %>%
  summarise(date = first(`date ringing`), age = age) %>%
  bind_cols(latitude = 52.51, longitude = 8.33, 
            site = "Dummer", country = "Germany")

gdu <- bind_rows(
  gdu[,c("metal", "date",  "age", "latitude", "longitude", "site", "country")],
  ringrecs
) %>% arrange(
  metal, date
)

##
gdu$combo <- NA
gdu$region <- NA
gdu$county <- NA
gdu$country <- NA
gdu$area <- NA
gdu$age <- NA

## observation type 
gdu$obstype <- NA

## scheme name
gdu$scheme <- "germ_JMelter"

## indicate these are color-ring data
gdu$obssource <- rep("cr")

##
gdu %<>% mutate(
  metal = as.character(metal),
  latitude = as.numeric(latitude),
  longitude = as.numeric(longitude)
)

## rmv obs w/out coords 
gdu %<>% filter(!is.na(latitude))
gdu %<>% filter(!is.na(longitude))

## remove few obs w/ lat/lon reversal for a few obs
gdu <- filter(gdu, !site %in% c("Hortobagy, Hajdu-Bihar", "Tjerkwerd, Hemdijk, NL", "Dingdener Heide", "Parque Natural de El Hondo, E"))

## get important columns for combining
gdu_clean <- gdu %>% dplyr::select(
  date, combo, metal, obstype, age, site, latitude, longitude, area,
  county, region, country, scheme, obssource)

## map it
sf::st_as_sf(gdu_clean,
             coords = c("longitude", "latitude"),
             crs = 4326, agr = "constant") %>%
  mapview::mapview()

## SAVE -----------------------------------------------------------------------

saveRDS(gdu_clean, "data/analysis/ringing/germany_JMelter_clean.rds")
