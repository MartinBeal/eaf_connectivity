## Alves team tracking data 

pacman::p_load(dplyr, sf, mapview, stringr, data.table, magrittr, amt)

## combine individuals stored w/ Piersma team data -----------------------------
folder <- "data/tracking/team_alves/islandica/"

list.files(folder, full.names = T)

alist <- lapply(list.files(folder, full.names = T), function(x){
  xx <- read.csv(x)
  xx <- dplyr::select(xx,
                      individual.local.identifier, timestamp, location.lat, location.long, 
                      sensor.type, argos.lc, study.name)
})

isl1 <- do.call(rbind, alist)

isl1$ageatrelease <- "adult" ## NEED TO CONFIRM THIS

isl1$study_site <- ifelse(
  isl1$individual.local.identifier %in% c("Kaldadarnes", "Vorsaber"), "Iceland",
  "NE")

## data from WaderTrack2020 ---------------------------------------------------
isl2  <- read.csv("data/tracking/team_alves/WaderTrack2020.csv")
meta <- read.csv("data/tracking/team_alves/WaderTrack2020-reference-data.csv")
meta <- subset(meta, animal.taxon == "Limosa limosa")

isl2 <- subset(isl2, individual.taxon.canonical.name == "Limosa limosa")

## satellite count filter
isl2 <- subset(isl2, gps.satellite.count > 3 | is.na(isl2$gps.satellite.count))

## study site 
isl2$study_site <- "PT"

## remove three indivs which only have local moves data in Tejo

isl2 %<>% 
  filter(!individual.local.identifier %in% c("J021767", "J021770", "J021730"))

## combine PTT and GPS datasets -----------------------------------------------
isl2$argos.lc <- "G"
isl2 %<>% dplyr::select(
  individual.local.identifier, timestamp, location.lat, location.long, 
  sensor.type, argos.lc, study.name, gps.hdop)

isl2$ageatrelease <- "adult" ## NEED TO CONFIRM THIS

isl <- bind_rows(isl1, isl2) 

isl$animal_sex <- "unknown"
isl$animal_ring_id <- "unknown" ## neither datasets include this...

## remove coord NAs
isl %<>% filter(!is.na(location.lat))

## remove Z locations (location process failure)
isl <- subset(isl, argos.lc != "Z")

## map it
# isl %>%
#    sf::st_as_sf(coords = c("location.long", "location.lat"), 
#                       crs = 4326, agr = "constant") %>% 
#    mapview()

## standardize for joining w/ other tracking datasets
isl <- isl %>% 
  mutate(
    timestamp = fasttime::fastPOSIXct(timestamp, tz = "UTC"),
    date      = fasttime::fastDate(timestamp),
  )

names(isl) <- gsub("\\.", "_", names(isl))
isl %<>% dplyr::rename(local_identifier=individual_local_identifier)


## order obs by id and time
isl <- data.table(isl[order(isl$local_identifier, isl$timestamp), ])

## SAVE ## -------------------------------------
saveRDS(isl, "data/analysis/tracking/Alves_trx_clean.rds")


## summarise sampling rate across devices -------------------------------------

one <- isl
# site  <- one$site_name[1]

one %<>% dplyr::select(-animal_sex, -date, -study_name)

## convert to amt 'tracks_xyt' format ##
tracks_amt <- one %>% 
  make_track(.x=location_long, .y=location_lat, .t=timestamp, 
             id = c(local_identifier, sensor_type), 
             crs = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_d"), 
             all_cols = T)

rate_summ <- tracks_amt %>%
  summarize_sampling_rate_many(
    cols = c("local_identifier", "sensor_type"), 
    time_unit = "hour")

ts_summ <- rate_summ %>% group_by(sensor_type) %>% 
  summarise(
    n_ids    = n_distinct(local_identifier),
    mn_md_ts = mean(na.omit(median)),
    mn_mn_ts = mean(mean),
    max_mn_tx = max(mean),
    min_mn_tx = min(mean)
  )

# ts_id_summ <- tracks_amt %>% group_by(season, id) %>% summarise(
#   md_ts = median(na.omit(ts))
# ) %>% ungroup() %>% summarise(
#   mn_md_id_ts = mean(md_ts)*60,
#   mx_md_id_ts = max(md_ts)*60
# )

