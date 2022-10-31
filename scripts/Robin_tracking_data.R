## process French btgo tracking data

## load packages
pacman::p_load(dplyr, magrittr, lubridate, stringr, mapview)

raw <- read.csv('data/tracking/team_robin/LIMITRACK [ID PROG 366].csv')

# raw %>% slice(1:10000) %>%
#   sf::st_as_sf(coords = c("location.long", "location.lat"),
#                      crs = 4326, agr = "constant") %>%
#   mapview()
# 
# raw %>% filter(tag.local.identifier %in% c("PAF25")) %>%
#   sf::st_as_sf(coords = c("location.long", "location.lat"),
#                crs = 4326, agr = "constant") %>%
#   mapview::mapview()

## w/ breeding migration
# "BQN01", "ISL13"
## uncertain
# MOZ16

raw %<>% group_by(tag.local.identifier) %>% 
  mutate(
    timestamp = fasttime::fastPOSIXct(timestamp, tz = "UTC")) 

# raw %>% filter(tag.local.identifier %in% "PAF25") %>% ggplot() +
#   geom_histogram(aes(x = (timestamp)))
# 
# raw %>% 
#   summarise(
#     start = first(timestamp),
#     end   = last(timestamp),
#     tdiff = difftime(end, start, units = 'days')
#   )

## Filter out birds which only had local movements (see Jourdan et al. 2022)
raw %<>% filter(tag.local.identifier %in% c("BQN01", "ISL13"))

raw %<>% filter(location.long < 180)

raw %<>% dplyr::select(
  individual.local.identifier, tag.local.identifier, location.long, 
  location.lat, timestamp, sensor.type
)


## get metal ring numbers (all numeric) from comments (most missing)

raw$animal_ring_id <- str_remove(
  raw$individual.local.identifier, pattern = fixed("["))
raw$animal_ring_id <- str_remove(
  raw$animal_ring_id, pattern = fixed("]"))
raw$animal_ring_id<- str_remove(
  raw$animal_ring_id, pattern = fixed("FRP_"))


raw_f <- raw %>% rename(
  location_long = location.long, location_lat = location.lat,
  local_identifier = tag.local.identifier, 
) %>% mutate(
  ageatrelease = "adult",
  gps.hdop = NA, argos_lc = "g", study_site = "Pertuis Charentais",
  animal_sex = NA, study_name = "Atl_France_islandica", sex = NA
) %>% dplyr::select(-individual.local.identifier)

# opd <- raw_f %>% group_by(local_identifier, date) %>% sample_n(1)
# opd <- raw_f %>% group_by(local_identifier, date) %>% summarise(timestamp=first(timestamp))


# tictoc::tic()
# opd <- raw_f[, .SD[sample(.N, 1)], by = list(local_identifier, date)]
# tictoc::toc()


## order obs by id and time
raw_f <- data.table::data.table(raw_f[order(raw_f$local_identifier, raw_f$timestamp), ])

## SAVE -----------------------------------------------------------------------
saveRDS(raw_f, "data/analysis/tracking/French_trx_clean.rds")


## summarise sampling rate across devices -------------------------------------

one <- raw_f
# site  <- one$site_name[1]

one %<>% dplyr::select(-animal_sex, -study_name)


## convert to amt 'tracks_xyt' format ##
tracks_amt <- one %>% 
  amt::make_track(.x=location_long, .y=location_lat, .t=timestamp, 
             id = c(local_identifier, sensor_type), 
             crs = sp::CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_d"), 
             all_cols = T)

rate_summ <- tracks_amt %>%
  amt::summarize_sampling_rate_many(
    cols = c("local_identifier"), 
    time_unit = "hour")


# ts_id_summ <- tracks_amt %>% group_by(season, id) %>% summarise(
#   md_ts = median(na.omit(ts))
# ) %>% ungroup() %>% summarise(
#   mn_md_id_ts = mean(md_ts)*60,
#   mx_md_id_ts = max(md_ts)*60
# )

