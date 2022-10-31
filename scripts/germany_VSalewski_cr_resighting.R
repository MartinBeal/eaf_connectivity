### Load German Schleswig-Holstein capture, recapture and resighting database - clean

## load packages
pacman::p_load(dplyr, magrittr, lubridate, stringr)

gsh <- readxl::read_xlsx('data/color/VSalewski_germany/BTG-Data_MOIN_2022-10-11.xlsx')

## translate column names 
gsh %<>% dplyr::rename(
  date = Datum, combo = `Colour rings`, metal = Ring, latitude = lat, 
  longitude = long, age = Age, site = Site, county = Province,
  country = Country
)

##
gsh$region <- NA
gsh$area <- NA

## rmv obs w/out coords 
gsh %<>% filter(!is.na(latitude))

## optional filter to remove records from breeding birds/probable breeders
# gsh %<>% filter(Status %in% c("non-breeding area", "first year"))

## observation type 
gsh$obstype <- NA

## scheme name
gsh$scheme <- "germ_VSalewski"

## indicate these are color-ring data
gsh$obssource <- rep("cr")

##
gsh %<>% mutate(
  metal = as.character(metal),
  age = as.character(age)
)

## get important columns for combining
gsh_clean <- gsh %>% dplyr::select(
  date, combo, metal, obstype, age, site, latitude, longitude, area,
  county, region, country, scheme, obssource)

gsh_clean <- filter(gsh_clean, longitude < 30)

## map it
# sf::st_as_sf(gsh_clean,
#              coords = c("longitude", "latitude"),
#              crs = 4326, agr = "constant") %>%
#   mapview::mapview()

## SAVE -----------------------------------------------------------------------

saveRDS(gsh_clean, "data/analysis/ringing/germany_VSalewski_clean.rds")
