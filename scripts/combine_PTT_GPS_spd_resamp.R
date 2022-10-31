## Combine GPS/PTT tracking data, speed filter, re-sample ---------------------

pacman::p_load(trip, dplyr, sf, amt, mapview)

dg  <- readRDS("data/analysis/tracking/Dutch_German_trx_clean.rds")
isl <- readRDS("data/analysis/tracking/Alves_trx_clean.rds")
ext <- readRDS("data/analysis/tracking/Gutierrez_trx_clean.rds")
fra <- readRDS("data/analysis/tracking/French_trx_clean.rds")

## combine datasets 

trx <- bind_rows(dg[, -"jdate"], isl, ext, fra)

trx %<>% dplyr::rename(
  ID = local_identifier, latitude = location_lat, longitude = location_long)

## Apply speed-angle filter ---------------------------------------------------

## convert to 'trip' object
trx %<>% dplyr::select(longitude, latitude, timestamp, ID, everything()) %>% 
  mutate(timestamp = format(timestamp))
tr <- trip::trip(trx, c("timestamp", "ID")) ## also acts as duplicate filter

n_ids <- n_distinct(tr$ID)
alist <- list()

for(i in seq_len(n_ids)){
  print(i)
  one <- subset(tr, ID == unique(tr$ID)[i])
  
  # one <- trip::trip(one, c("timestamp", "ID")) ## also acts as duplicate filter
    
  one$mcfilter <- trip::speedfilter(one, max.speed = 80)    # Create a filter of a track for "bad" points implying a speed of motion that is unrealistic.
  # print(one@proj4string)
  alist[[i]] <- one
}

trx <- do.call(rbind, alist)

### visualize
# mapview(head(trx, 10000), zcol = "mcfilter")
# mapview(tail(trx, 10000), zcol = "mcfilter")
# mapview(subset(trx, study_name == "Black-tailed Godwits Extremadura (Lotek)"), 
#         zcol = "mcfilter")

## filter -----------
trx_f <- as.data.frame(subset(trx, mcfilter == TRUE))


## summarise sampling rate across devices -------------------------------------

one <- trx_f
# site  <- one$site_name[1]

one %<>% dplyr::select(-animal_sex, -date, -study_name)


## convert to amt 'tracks_xyt' format ##
tracks_amt <- one %>% 
  make_track(.x=longitude, .y=latitude, .t=timestamp, 
             id = c(ID), 
             crs = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_d"), 
             all_cols = T)

rate_summ <- tracks_amt %>%
  summarize_sampling_rate_many(
    cols = c("ID"), 
    time_unit = "hour")

ts_summ <- rate_summ %>%
  summarise(
    n_ids    = n_distinct(ID),
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


## re-sample to 12h/24h -------------------------------------------------------

trk1 <- tracks_amt %>% nest(data = -"ID")

trk2 <- trk1 %>% 
  mutate(steps = map(data, function(x) 
    x %>% track_resample(rate = hours(12), tolerance = hours(6)))
  )

trk2 <- trk2 %>% select(ID, steps) %>% unnest(cols = steps)
trk2

## re-name things 

trk2 %<>% dplyr::rename(
  id = ID, longitude = x_, latitude = y_, timestamp = t_, 
  device = sensor_type, age = ageatrelease) %>% dplyr::select(-burst_) %>% 
  as_tibble()

## SAVE ## --------------------------------------------------------------------

saveRDS(trk2, "data/analysis/tracking/PTT_GPS_mconn_12h.rds")
