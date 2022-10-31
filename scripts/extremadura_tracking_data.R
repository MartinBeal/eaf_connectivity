## Tracking data from Extremadure

## data paid for and collected by Extremadura team
raw <- read.csv("data/tracking/team_gutierrez/Black-tailed Godwits Extremadura (Lotek).csv")
meta <- read.csv("data/tracking/team_gutierrez/Black-tailed Godwits Extremadura (Lotek)-reference-data.csv") 

## from movebank:
# Study Reference Location	
# Longitude	-5.975
# Latitude	39.038

## All individuals have both PTT (only lc 1,2,3) and GPS data 
raw$argos.lc <- ifelse(is.na(raw$argos.lc), "G", raw$argos.lc)

## combine PTT and GPS datasets -----------------------------------------------
raw %<>% dplyr::select(
  individual.local.identifier, timestamp, location.lat, location.long, 
  sensor.type, argos.lc, study.name)

raw$gps.hdop <- NA # not provided

raw$ageatrelease <- "adult" ## NEED TO CONFIRM THIS

raw$animal_sex <- "f"

raw %<>% left_join(meta[c("animal.ring.id", "animal.id")], 
                   by=c("individual.local.identifier"="animal.id"))

## remove coord NAs
raw %<>% filter(!is.na(location.lat))

## broad bounding box filter
raw %<>% filter(location.long > -30 & location.long < 26)
raw %<>% filter(location.lat > 7    & location.lat < 58)

## map it
# raw %>%
#    sf::st_as_sf(coords = c("location.long", "location.lat"),
#                       crs = 4326, agr = "constant") %>%
#    mapview()

## standardize for joining w/ other tracking datasets
raw <- raw %>% 
  mutate(
    timestamp = fasttime::fastPOSIXct(timestamp, tz = "UTC"),
    date      = fasttime::fastDate(timestamp),
  )

names(raw) <- gsub("\\.", "_", names(raw))
raw %<>% rename(local_identifier=individual_local_identifier)

## remove indivs w/ few locations (<=3)
keepids <- raw %>% group_by(local_identifier) %>% summarise(n_pnts = n()) %>% 
  filter(n_pnts > 3)
raw %<>% filter(local_identifier %in% keepids$local_identifier)

## order obs by id and time
raw <- data.table(raw[order(raw$local_identifier, raw$timestamp), ])

## remove data w/ bad time stamps (not from 2022)
raw <- subset(raw, year(raw$timestamp) == 2022)

raw$study_site <- "Iberia"

## SAVE -----------------------------------------------------------------------
saveRDS(raw, "data/analysis/tracking/Gutierrez_trx_clean.rds")


## summarise sampling rate across devices -------------------------------------
one <- raw
# site  <- one$site_name[1]

one %<>% dplyr::select(-animal_sex, -date, -study_name)

## convert to amt 'tracks_xyt' format ##
tracks_amt <- one %>% 
  make_track(.x=location_long, .y=location_lat, .t=timestamp, 
             id = c(local_identifier), 
             crs = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_d"), 
             all_cols = T)

rate_summ <- tracks_amt %>%
  summarize_sampling_rate_many(
    cols = c("local_identifier"), 
    time_unit = "hour")

ts_summ <- rate_summ %>% 
  summarise(
    n_ids    = n_distinct(local_identifier),
    mn_md_ts = mean(na.omit(median)),
    mn_mn_ts = mean(mean),
    max_mn_tx = max(mean),
    min_mn_tx = min(mean)
  )
ts_summ

# ts_id_summ <- tracks_amt %>% group_by(season, id) %>% summarise(
#   md_ts = median(na.omit(ts))
# ) %>% ungroup() %>% summarise(
#   mn_md_id_ts = mean(md_ts)*60,
#   mx_md_id_ts = max(md_ts)*60
# )



## re-sample to 12h/24h -------------------------------------------------------

trk1 <- tracks_amt %>% nest(data = -"local_identifier")

trk2 <- trk1 %>% 
  mutate(steps = map(data, function(x) 
    x %>% track_resample(rate = hours(12), tolerance = hours(6)))
  )

trk2 <- trk2 %>% select(local_identifier, steps) %>% unnest(cols = steps)
trk2

## re-name things 

trk2 %<>% dplyr::rename(
  id = local_identifier, longitude = x_, latitude = y_, Datetime = t_, 
  device = sensor_type, age = ageatrelease) %>% dplyr::select(-burst_) %>% 
  as_tibble()

## map it
trk2 %>%
  sf::st_as_sf(coords = c("longitude", "latitude"),
               crs = 4326, agr = "constant") %>%
  mapview()

## SAVE ##