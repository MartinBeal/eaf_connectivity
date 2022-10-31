## Streamline dataset by removing observations of the same individual closely
# spaced in time and space

pacman::p_load(data.table, dplyr, magrittr)

alldat <- readRDS("data/analysis/ringing/comb_euring_cring.rds")


## get time difference between obs per bird -----------------------------------

tictoc::tic()
# creating an ordered data.table
alldat <- data.table(alldat[order(alldat$bird_id, alldat$date), ])

# super fast:
# source: (https://stackoverflow.com/questions/32999460/how-to-calculate-time-difference-with-previous-row-of-a-data-frame-by-group)
alldat[ , tdiff := c(0, diff(as.numeric(date))), by=bird_id]
tictoc::toc()

## use 123456789 to replace NAs (an annoyance in ifelse statements to come)
alldat$tdiff <- ifelse(as.numeric(alldat$tdiff) == 0, 123476789, alldat$tdiff)


## calculate distance between subsequent locations w/in individuals -----------

# need to remove locations w/out lat long first
alldat <- subset(alldat, !is.na(latitude))
alldat <- subset(alldat, !is.na(longitude))

## again remove birds w/ only one sighting
x <- alldat %>% group_by(bird_id) %>% summarise(nobs = n()) %>% filter(nobs == 1)
alldat <- filter(alldat, !bird_id %in% x$bird_id)

## assure time-bird order
setorder(alldat, bird_id, date)

## calculate distance (data.table way, fast)
alldat[, dist := c(NA,
  round(geodist::geodist(
    cbind(longitude, latitude),
    measure='haversine', sequential = TRUE)/1000,1)),
  by = bird_id]

## use 123456789 to replace NAs (an annoyance in ifelse statements to come)
alldat$dist <- ifelse(is.na(alldat$dist), 123456789, alldat$dist)


## Identify when consecutive obs are from same 'site' -------------------------

## get duplicated consecutive obs, per bird (data.table way)
alldat[, samecnt := data.table::rleid(site),
       by = bird_id]
alldat$samecnt <- ifelse(is.na(alldat$site), NA, alldat$samecnt) # NA where site NA

## get logical vector identifying first of a string of obs at a site 
alldat %<>%
  group_by(temp_id = paste(bird_id, samecnt)) %>%
  mutate(first = as.integer(row_number() == 1L)) %>%
  ungroup() %>% 
  mutate(
    isfirst = ifelse(is.na(alldat$site), NA, first) # NA where site NA
  )

## Filter ---------------------------------------------------------------------
## #1. identify obs w/in a certain time window (e.g. a week)
## #2. remove recent reobs at the same site OR less than a distance thresh 
bdist <- 20 # km   - buffer distance
ttime <- 7 # days - time window threshold

## ID obs w/in the reobs time limit
alldat$shortreobs <- ifelse(
  alldat$tdiff < ttime, TRUE, FALSE
)

## among short obs, which are close to last obs?
alldat$tormv <- ifelse(
  alldat$shortreobs == TRUE & (
    (is.na(alldat$site) & alldat$dist < bdist) | 
      (!is.na(alldat$site) & alldat$isfirst == 0)
  ), TRUE, FALSE
)

## do some manual checks
# View(dplyr::select(alldat, bird_id, shortreobs, site, isfirst, tdiff, dist, tormv))

alldat_f2 <- subset(alldat, tormv == FALSE)

## again remove birds w/ only one sighting
xz <- alldat_f2 %>% group_by(bird_id) %>% summarise(nobs = n()) %>% filter(nobs == 1)
alldat_f2 <- subset(alldat_f2, !bird_id %in% xz$bird_id)

## clean it up
alldat_f2 %<>% 
  dplyr::select(bird_id, combo, metal, date, longitude, latitude, obstype, 
                obssource, site, area, country, scheme_country)

## SAVE -----------------------------------------------------------------------

saveRDS(alldat_f2,
        paste0("data/analysis/ringing/comb_euring_cring_no", ttime,
               "dayreobs.rds")
        )
