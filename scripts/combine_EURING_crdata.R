### Combine color-ringing and EURING datasets ---------------------------------

pacman::p_unload(plyr)
pacman::p_load(dplyr, magrittr, stringr)

# crdat <- readRDS("data/analysis/ringing/cr_tag_ice_wash_merge.rds") # from cr-schemes
crdat <- readRDS("data/analysis/ringing/cr_merge.rds") # from cr-schemes

ll <- readRDS("data/analysis/ringing/EURING_BTGO.rds") # limosa limosa

## Convert cr-schemes into countries
crschemes <- read.csv("data/color/scheme_countries.csv")
crdat     <- left_join(crdat, crschemes)

## Identifying shared records btwn EURING and scheme-leader datasets

## how many metal ids are shared?
euids <- unique(ll$metal)
crids <- unique(crdat$metal)

sum(euids %in% crids)
sum(crids %in% euids)

# get ids w/ data in both datasets 
ids <- euids[which(euids %in% crids)]

crdat$date <- as.Date(crdat$date)

crdat %<>% as_tibble()

# ll_f    <- anti_join(ll, crdat, by=c("bird_id", "metal", "date")) # rmv shared obs from euring data
# shared  <- semi_join(crdat[,c(1,2,4)], ll[,c(1,2,4)], by=c("bird_id", "metal", "date")) #id them
ll_f    <- anti_join(ll, crdat, by=c("metal", "date")) # rmv shared obs from euring data
shared  <- semi_join(crdat[,c(1,3)], ll[,c(2,4)], by=c("metal", "date")) #id them

alldat  <- full_join( # join filtered euring data w/ all crdata
  crdat, ll_f, 
  by=c("combo", "metal", "date", "longitude", "latitude", "obstype", "obssource", "scheme_country")) # combine
# alldat  <- full_join(crdat, ll, by=c("metal", "date")) # combine

## clean it up
alldat %<>% dplyr::select(
  "combo", "metal", "date", "longitude", "latitude", "obstype", 
  "obssource", "site", "area", "county", "country", "region", "scheme_country")

## create universal bird_id w/ either metal or combo or both

alldat$bird_id <- ifelse( # combine cr and metal code, use _
  is.na(alldat$combo),    # to signify where one or other is missing
  alldat$metal,
  ifelse(
    is.na(alldat$metal),
    alldat$combo,
    paste(alldat$combo, alldat$metal, sep="_")
  )
)

## remove all birds w/ only one obs (assuming these are true non-resighted brds)
xx <- alldat %>% group_by(bird_id) %>% summarise(nobs = n()) %>% filter(nobs == 1)
xxx <- subset(alldat, bird_id %in% xx$bird_id)
table(xxx$scheme_country) # which schemes/countries have ringing only data?

alldat <- subset(alldat, !bird_id %in% xx$bird_id)

## remove locations with no coordinates ---------------------------------------

alldat %<>% filter(!is.na(latitude))
alldat %<>% filter(!is.na(longitude))

## Save -----------------------------------------------------------------------

saveRDS(alldat, "data/analysis/ringing/comb_euring_cring.rds")
