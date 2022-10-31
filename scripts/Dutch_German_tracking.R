
pacman::p_load(dplyr, sf, mapview, stringr, data.table, magrittr, amt)


## Dutch team tracking data 
raw  <- read.csv("data/tracking/team_piersma/allLocations.csv")
meta <- read.csv("data/tracking/team_piersma/combinedReferenceData.csv")


raw$argos_lc <- ifelse(raw$sensor_type == "GPS", "G", raw$argos_lc)

## Remove hand-raised birds 

raw <- raw[!str_detect(raw$study_name, pattern="hand-raised"), ]

gps <- subset(raw, sensor_type == "GPS")
ptt <- subset(raw, sensor_type != "GPS")

# gps %>%
#   sf::st_as_sf(coords = c("location_long", "location_lat"), 
#                      crs = 4326, agr = "constant") %>% 
#   mapview()
# 
# ptt %>% 
#   filter(local_identifier %in% unique(ptt$local_identifier)[1:10]) %>%
#   sf::st_as_sf(coords = c("location_long", "location_lat"), 
#                crs = 4326, agr = "constant") %>% 
#   mapview()

iber <- subset(raw, study_site == "Iberia")

ifir <- iber %>% group_by()

## Identify birds which were definitely juveniles at release 

raw$ageatrelease <- ifelse(
  str_detect(raw$study_name, pattern="juveniles") | 
                          str_detect(raw$study_name, pattern="chicks"),
  "juvenile", "adult"
)

raw %>% group_by(ageatrelease) %>% summarise(n_id = n_distinct(local_identifier))

## get metal ring numbers (all numeric) from comments (most missing)

raw$metal_ring_id <- ifelse(
  str_detect(raw$comments, "^[:digit:]+$"), raw$comments, NA)

## when ring number missing from animal_ring_id bring (numeric code) 
# from local_identifier

raw$animal_ring_id <- ifelse(
  is.na(raw$animal_ring_id) & str_detect(raw$local_identifier, "^[:digit:]+$"), 
  raw$local_identifier, raw$animal_ring_id)

## get 1/2 points per day per individual

raw <- raw %>% 
  mutate(
    timestamp = fasttime::fastPOSIXct(raw$timestamp, tz = "UTC"),
    date      = fasttime::fastDate(raw$timestamp),
    jdate     = round(julian.Date(date),0)
  )

raw_f <- dplyr::select(
  raw, local_identifier, study_name, timestamp, sensor_type,jdate, date, 
  ageatrelease, location_lat, location_long, animal_ring_id, animal_sex, 
  study_site, argos_lc, gps_hdop
  )

# opd <- raw_f %>% group_by(local_identifier, date) %>% sample_n(1)
# opd <- raw_f %>% group_by(local_identifier, date) %>% summarise(timestamp=first(timestamp))


# tictoc::tic()
# opd <- raw_f[, .SD[sample(.N, 1)], by = list(local_identifier, date)]
# tictoc::toc()

## remove indivs w/ few locations (<10)
keepids <- raw_f %>% group_by(local_identifier) %>% summarise(n_pnts = n()) %>% 
  filter(n_pnts > 10)

raw_f %<>% filter(local_identifier %in% keepids$local_identifier)

## order obs by id and time
raw_f <- data.table(raw_f[order(raw_f$local_identifier, raw_f$timestamp), ])

## SAVE -----------------------------------------------------------------------
saveRDS(raw_f, "data/analysis/tracking/Dutch_German_trx_clean.rds")


## summarise sampling rate across devices -------------------------------------

one <- raw_f
# site  <- one$site_name[1]

one %<>% dplyr::select(-jdate, -animal_sex, -date, -study_name)


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

