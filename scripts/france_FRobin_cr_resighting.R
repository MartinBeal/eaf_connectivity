### Load French capture, recapture and resighting database - clean, combine

## load packages
pacman::p_load(dplyr, magrittr, lubridate, stringr)

fren <- read.delim('data/color/FRobin_france/BaseLimicoles_LIMLIM_ExpBase_20221007_024112.csv', sep=";")

## Limosa limosa
fren <- subset(fren, espece == "LIMLIM")

## non-breeding ringing in Moeze Oleron and limosa subsp chicks and adults
fren <- subset(fren, id_prog_bag %in% c("366", "600", "SENE")) # maybe "SENE" too?

## keep only validated obs
fren %<>% filter(validation == "oui")

## keep only records under French scheme
fren %<>% filter(centre == "FRP")

## translate column names 
fren %<>% dplyr::rename(
  date = date_obs, combo = coul_pat_combi, metal = bague, 
  site = lieudit, area = localite, county = dept, country = pays
)

##
fren$region <- NA

##
fren %<>% filter(latitude != 0)

## date to datetime
fren$date <- as.Date(fren$date, format="%d/%m/%Y")

## scheme name
fren$scheme <- "fran_FRobin"

## indicate these are color-ring data
fren$obssource <- rep("cr")

## observation type 
fren$obstype <- ifelse(fren$action == "B", "N", 
                       ifelse(fren$action == "C", "S", fren$action))

## birds missing metal rings
fren$metal <- ifelse(fren$metal == "", NA, fren$metal)

## birds only ringed, never re-seen
x <- fren %>% group_by(metal) %>% summarise(n = n()) %>% filter(n==1)
fren <- filter(fren, !metal %in% x$metal)

## get important columns for combining
fren_clean <- fren %>% dplyr::select(
  date, combo, metal, obstype, age, site, latitude, longitude, area,
  county, region, country, scheme, obssource)

## map it
fren_sf <- sf::st_as_sf(fren_clean,
             coords = c("longitude", "latitude"),
             crs = 4326, agr = "constant") 
fren_sf %>%
  mapview::mapview()

## custom removal of points in clearly wrong place (e.g. over sea)
fren_clean <- fren_clean[-which(fren_clean$site == "wwt wetland centre" & fren_clean$area == "SLIMBRIDGE")
, ]
fren_clean <- fren_clean[-which(fren_clean$site == "wildfowl and wetlands trust" & fren_clean$area == "LLANELLI"), ]

fren_clean$latitude <- ifelse(fren_clean$site == "parque natural de l'albufera" & fren_clean$area == "EL PALMAR", 39.312200, fren_clean$latitude)
fren_clean$longitude <- ifelse(fren_clean$site == "parque natural de l'albufera" & fren_clean$area == "EL PALMAR", -0.326385, fren_clean$longitude)

fren_clean <- fren_clean[-c(1338, 1345, 1380:1396, 1574, 1576, 14036, 15465, 16040, 17920, 21089), ]
fren_clean <- fren_clean[-c(17889, 17890), ]

## SAVE -----------------------------------------------------------------------

saveRDS(fren_clean, "data/analysis/ringing/france_FRobin_clean.rds")
