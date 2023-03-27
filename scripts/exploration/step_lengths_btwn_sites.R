## get step length (btwn sites) distribution from tracking data

pacman::p_load(amt, data.table)


## which network to create
# season <- "all"
# season <- "spring"
season <- "fall"

## Load combined data, with both IBA overlay info and hexgrid outsites -------
### SAVE 
alldat <- readRDS(
  paste0("data/analysis/combined/alldatatypes_ibas_outsites_all.rds"))

trax <- filter(alldat, datatype == "trax")

if(season == "all"){ ## year-round
  trax <- trax ## all individuals
} else if(season == "spring"){
  trax <- subset(trax, lubridate::month(trax$timestamp) %in% c(1:6))
} else if(season == "fall"){
  doy <- lubridate::yday(trax$timestamp) # June 23/24: 175
  trax <- subset(
    trax, 
    month(trax$timestamp) %in% c(1, 7:12) | doy %in% c(175:181)
  )
}


## remove birds w/ only one sighting (to avoid unconnected nodes) -------
xz <- trax %>% group_by(bird_id) %>% 
  summarise(nobs = n()) %>% filter(nobs == 1)
trax <- subset(trax, !bird_id %in% xz$bird_id)

## remove birds only seen (multiple times) at same site ----------------------- 
nsites <- trax %>% group_by(bird_id) %>% 
  summarise(nsites = n_distinct(SitRecID))

## % of birds w/ relocs at one site only
sum(nsites$nsites == 1) / n_distinct(trax$bird_id) * 100

xz2 <- filter(nsites, nsites == 1)
trax <- subset(trax, !bird_id %in% xz2$bird_id)


### one point per visit to/period at a site
trax$SitRecID2 <- trax$SitRecID
trax_f <- as.data.table(trax)[, .SD[1], by = list(bird_id, rleid(SitRecID2))]

## convert to amt 'tracks_xyt' format ##
tracks_amt <- trax_f %>% 
  make_track(.x=longitude, .y=latitude, .t=timestamp, 
             id = c("bird_id"), 
             crs = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_d"), 
             all_cols = T) %>% 
  arrange(bird_id, t_)

trax_f$step_length <- step_lengths(tracks_amt, lonlat=T, append_last = T) / 1000

trckstep <- steps(tracks_amt, lonlat=T, append_last = T)

trckstep <- split(tracks_amt, tracks_amt$bird_id) %>% 
  # map(steps(lonlat=T, append_last = T)) 
  map(function(x) {
    # print(x$bird_id[1])
    xx <- steps(x, lonlat=T, append_last = T)
    xx$bird_id <- x$bird_id[1]
    return(xx)
    }) %>% 
  bind_rows() %>% 
  mutate(
    sl_ = sl_/1000, # km
    step_time = round(as.numeric(dt_)/60/60, 2), # hours
    direction_deg = REdaS::rad2deg(direction_p),
  ) %>% rename(
    direction_rad = direction_p
  )

hist(trckstep$step_time)

hist(subset(trckstep, step_time <= 168)$sl_, 50)

hist(subset(trckstep, step_time <= 336)$sl_, 50)
hist(subset(trckstep, step_time <= 280)$sl_, 50)
hist(subset(trckstep, sl_ > 1500)$step_time, 50)

hist(subset(trckstep, sl_/1000 > 1500)$step_time, 50)

## direction
hist(trckstep$direction_rad)

## split by latitude grouping 

trckstep$lat_group <- as.numeric(cut(trckstep$y1_, 4)) # equal-length
trckstep$lon_group <- as.numeric(cut(trckstep$x1_, 4)) # equal-length

trckstep$latlon_group <- paste(trckstep$lat_group, trckstep$lon_group)

ggplot() + 
  # geom_histogram(data = trckstep, aes(x=direction_rad)) +
  # geom_histogram(data = trckstep, aes(x=direction_deg)) +
  geom_histogram(data = trckstep, aes(x=sl_)) +
  facet_wrap(~latlon_group)

trckstep %>% group_by(lat_group) %>% 
  summarise(
    med_displ = median(sl_)
  )

##

# trckstep %<>% filter(sl_ > 10)

hist(trckstep$sl_, 50, main=NULL)

med_displ <- median(trckstep$sl_)
max_displ <- max(trckstep$sl_)


