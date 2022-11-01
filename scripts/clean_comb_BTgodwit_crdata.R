### Load all data and format consistently for analysis
## Edited by Martin Beal to include all resighting (not just in Tagus)

## load packages
pacman::p_load(dplyr, magrittr, plyr, lubridate)


#### renaming functions #### MB Add metal ring combo (to combing w/)
# caps_names <- function(x) {
#   names(x) <- c('combo', 'capdate', 'site', 'age')
#   return(x)
# }
# dat_names <- function(x) {
#   names(x) <- c('combo', 'date', 'site')
#   return(x)
# }
re_names <- function(x) {
  names(x) <- c('metal', 'combo', 'date', 'site', 'age', 'obstype', 'latitude', 
                'longitude', 'area', 'county', 'region', 'country')
  return(x)
}


# lookup table of locations, for renaming and connecting to coords
site_lookup <- read.csv('data/godwitdata_josh/site_lookup_airport.csv')
# site_lookup <- read.csv("data/godwitdata_josh/regodwitdata/lookup.csv")

rename_sites <- function(x) {
  mysite <- site_lookup$mysite[match(x, site_lookup$oldsite)]
}


#### TAGUS DATASET #### -------------------------------------------------------
# tag <- read.csv('data/godwitdata/PortugalBase06.21.csv')[,-1]#, fileEncoding = "UTF-16")
tag <- readxl::read_xlsx('data/godwitdata_josh/PortugalBase06.21.xlsx')[,-1]#, fileEncoding = "UTF-16")

## observer names (many repeated, and lists of multiple peeps)
sort(unique(tag$`Resight Observer`))

## "status"is New (ie, ringed for the first time), Recapture, Sighting or X=dead)
# formatting
tag$Date %<>% as.POSIXct(format='%d.%m.%y')
tag$doy <- as.numeric(format(tag$Date, '%j')) # day of year
names(tag)[7] <- 'combo'

### Fill out metal #s and identify when combo re-used after death -------------
# which cr-combos have multiple metal ring combos?
xx <- subset(tag, !is.na(tag$Metal))
reuse_ids <- xx %>% dplyr::group_by(combo) %>% 
  dplyr::summarise(n_metal=n_distinct(Metal)) %>% filter(n_metal > 1)

## how many obs have no info on type of observation?
sum(is.na(tag$Status))

tag <- subset(tag, !is.na(tag$Status)) # filter to only obs w/ obstype

n_crs <- n_distinct(tag$combo) # # of cr-combos
## loop thru cr-combos, filling out metal combos where missing
alist <- lapply(seq_len(n_crs), function(x){
  print(x)
  one <- subset(tag, combo == unique(tag$combo)[x])
  n_met <- n_distinct(na.omit(one$Metal))
  
  if(n_met == 1){
    one$Metal <- unique(na.omit(one$Metal))
  } else if(n_met > 1) {
    one <- arrange(one, one$Date)
    ids_met <- unique(na.omit(one$Metal))
    
    k <- 1
    for(i in seq_len(nrow(one))){
      rw <- one[i, ]
      id <- ids_met[k] # get metal combo of living bird (at that time)
      
      if(is.na(rw$Metal)) { one$Metal[i] <- id } # fill in NAs for metal ring
      
      if(rw$Status == "X") {k <- k + 1} # if bird died, move to next metal combo
      
    }
  }
  return(one)
})

rcmbn <- do.call(rbind, alist)
# check that the filling-in worked for specific cases
# View(subset(rcmbn, combo == reuse_ids$combo[1]))

tag <- rcmbn

## some records don't have a site
nosite <- which(is.na(tag$tbl_location.site))
tag$tbl_location.site[nosite] <- tag$Site...32[nosite]
tag$tbl_location.site %<>% factor

## use multiple location columns to infer which is correct and fill in NAs
tag$site <- ifelse(is.na(tag$Site...32), tag$tbl_location.site, tag$Site...32)
tag$region <- ifelse(is.na(tag$Region...28), tag$tbl_location.region, tag$Region...28)
tag$area <- tag$Area
tag$county <- tag$County
tag$country <- tag$`Resight country`

## all caps and sightings, clean
# location info at three scales (small to large): site, area, region
tag_red <- tag %>% 
  dplyr::select(Metal, combo, Date, site, Age, Status, n, w, area, county, 
                region, country) %>% 
  re_names

## add coordinates for sites in Tejo missing them
tag_red <- left_join(tag_red, site_lookup[,c("oldsite", "n", "w")], by=c("site" = "oldsite"))

tag_red$latitude  <- as.numeric(tag_red$latitude)
tag_red$longitude <- as.numeric(tag_red$longitude)
tag_red$latitude  <- ifelse(is.na(tag_red$latitude), tag_red$n, tag_red$latitude)
tag_red$longitude <- ifelse(is.na(tag_red$longitude), tag_red$w, tag_red$longitude)

tag_clean <- tag_red %>% 
  dplyr::select(metal, combo, date, site, age, obstype, latitude, longitude, 
                area, county, region, country)

tag_clean$scheme <- "tag_JAlves"


# re-captures
# caps_tag <- subset(tag, Status == 'N', select=c(combo, Date, tbl_location.site, Age)) %>% 
  # droplevels %>% caps_names

# re-sightings only, # 2000 onwards
# dat_tag <- subset(tag, Status=='S', select=c(combo, Date, tbl_location.site)) %>% 
  # droplevels %>% dat_names #%>% 
  #subset(date > as.POSIXct('2000-01-01') & date < as.POSIXct('2021-06-01'))



#### ICELAND DATASET #### -----------------------------------------------------
# ice <- read.csv('data/godwitdata/ice_godwit_tagus.30.06.2021.csv')
ice <- readxl::read_xlsx('data/godwitdata_josh/Iceland.30.03.2022.xlsx')#, fileEncoding = "UTF-16")
ice$Newdate %<>% as.POSIXct(format='%d/%m/%Y')

# format
ice <- dplyr::rename(
  ice, 
  "combo" = "Colour-ring combination", "Metal"="ringnumber")
ice$combo %<>% toupper() # uppercase, for consistency of combinations

## remove one individual w/ two single obs w/ different color combos
id <- ice %>% 
  dplyr::group_by(Metal) %>% 
  dplyr::summarise(n=n_distinct(combo)) %>% dplyr::filter(n>1)

ice %<>% dplyr::filter(Metal != id$Metal)

## where no site specified, indicate there are coordinates 
ice$site <- ifelse(is.na(ice$site), "coords", ice$site)

# extract captures - capture info not included :()
# caps_ice <- subset(ice, observartion=='Captured', select=c(combo, Newdate, site, age)) %>% 
# droplevels %>% caps_names

ice %<>% mutate(latitude = NA, longitude = NA)
ice_clean <- ice %>% 
  dplyr::select(Metal, combo, Newdate, site, age, observartion, n, w, area, 
                county, region, country) %>% 
  re_names

## standardize Status/obstype
ice_clean$obstype <- ifelse(
  ice_clean$obstype == "Captured", "N",
  ifelse(ice_clean$obstype == "Dead", "X", 
         ifelse(ice_clean$obstype == "Recaptured", "R", 
                ifelse(ice_clean$obstype %in% c("Sighted", "sighted"), "S",
                       ice_clean$obstype))))

ice_clean$scheme <- "ice_TGunn"


#### WWRG dataset #### -------------------------------------------------------
# wash <- read.csv('data/godwitdata/Wash CR Godwits.csv')
wash <- readxl::read_xlsx('data/godwitdata_josh/Wash CR Godwits.xlsx')

# formatting
wash <- dplyr::rename(wash, combo = "Colour-rings", Metal = `Ring No.`)

# combine year and month (w/ random day) for those obs missing Date row
wash$Date <- ifelse(
  is.na(wash$Date), 
  paste(wash$Year, wash$month...44, round(runif(1, 1, 30), 0), sep = "-"), 
  wash$Date)

wash$Date %<>% as.POSIXct(format='%d.%m.%y')

### Fill out metal #s and identify when combo re-used after death -------------
## remove few obs w/out cr-combo
wash <- subset(wash, !is.na(combo))

# which cr-combos have multiple metal ring combos? **NONE in this dataset
# xx <- subset(wash, !is.na(wash$Metal))
# reuse_ids <- xx %>% dplyr::group_by(combo) %>% 
#   dplyr::summarise(n_metal=n_distinct(Metal)) %>% 
#   dplyr::filter(n_metal > 1)

## how many obs have no info on type of observation: NONE?
sum(is.na(wash$Status))

## remove few obs w/out cr-combo
wash <- subset(wash, !is.na(combo))

n_crs <- n_distinct(wash$combo) # # of cr-combos
## loop thru cr-combos, filling out metal combos where missing
alist <- lapply(seq_len(n_crs), function(x){
  print(x)
  one <- subset(wash, combo == unique(wash$combo)[x])
  n_met <- n_distinct(na.omit(one$Metal))
  
  if(n_met == 1){
    one$Metal <- unique(na.omit(one$Metal))
  }
  return(one)
})

rcmbn <- do.call(rbind, alist)
# check that the filling-in worked for specific cases
# View(subset(rcmbn, combo == unique(rcmbn$combo[1])))

wash <- rcmbn

## some sightings have mismatched info...
# View(wash[which(wash$country != wash$`Resight country`),
#           c("country", "Resight country", "n", "w")])

wash$country <- ifelse(wash$country == "Morroco", "Morocco", wash$country)
wash$country <- ifelse(wash$country == "Channel Island", "Channel Islands", wash$country)
wash$`Resight country` <- ifelse(
  wash$country == "Northern Ireland" & wash$`Resight country` == "Ireland",
  "Northern Ireland", wash$`Resight country`)
wash$`Resight country` <- ifelse(
  wash$country == "Wales" & wash$`Resight country` == "England",
  "Wales", wash$`Resight country`)

wash[which(wash$country != wash$`Resight country`),
          c("country", "Resight country", "n", "w")]

## remove the leftovers (too wrong to know which is right!)
wash <- wash[-which(wash$country != wash$`Resight country`), ]

## fix latitude values erroneously above 99 deg
wash[which(wash$n > 99),"n"] <- as.numeric(substring((subset(wash, wash$n > 99)$n), 2))

## use multiple location columns to infer which is correct and fill in NAs
wash$site <- ifelse(is.na(wash$tbl_location.site), wash$Site...31, wash$tbl_location.site)

wash$area <- wash$Area

wash$region <- ifelse(
  wash$tbl_location.region == wash$Region...27, 
  wash$tbl_location.region,
  ifelse(wash$tbl_location.region == wash$Region...50, 
         wash$tbl_location.region, 
         wash$Region...27))

wash$county <- ifelse(is.na(wash$County), wash$tbl_location.county, wash$County)

wash$country <- wash$`Resight country`

## clean it up
wash_clean <- wash %>% 
  dplyr::select(Metal, combo, Date, site, Age, Status, n, w, area, county, 
                region, country) %>% 
  re_names

wash_clean$scheme <- "wash_JGill"

## separating captures and re-sightings 
# captures are those records that include ring number
# caps_wash <- subset(
#   wash, 
#   !is.na("Ring No."), 
#   select=c(combo, Date, tbl_location.site, Age)) %>% 
#   droplevels %>% 
#   caps_names

# sightings from the Tagus estuary
# dat_wash <- subset(
#   wash, tbl_location.area == 'Tagus estuary', 
#   select=c(combo, Date, tbl_location.site)) %>% 
#   droplevels %>% 
#   dat_names %>% 
#   subset(date > as.POSIXct('2000-01-01')  & date < as.POSIXct('2021-06-01'))


#### Combine datasets ####

# caps <- rbind(caps_tag, caps_wash) # NB no Ice
# caps$site %<>% rename_sites()

# dat <- rbind(dat_tag, dat_ice, dat_wash)
# dat$site %<>% trimws %>% rename_sites() ## trim whitespace

## how many obs are missing site info?
alist <- list(tag_clean, ice_clean, wash_clean)

lapply(alist, function(x) sum(is.na(x$site)))
lapply(alist, function(x) sum(is.na(x$area)))
lapply(alist, function(x) sum(is.na(x$region)))
lapply(alist, function(x) sum(is.na(x$county)))
lapply(alist, function(x) sum(is.na(x$country)))
lapply(alist, function(x) sum(is.na(x$latitude)))
lapply(alist, function(x) sum(is.na(x$date)))

## combine 
alldat <- rbind(tag_clean, ice_clean, wash_clean)

## remove obs w/out date info
alldat %<>% filter(!is.na(date))

## rename sites by lookup table (WHICH IS INCOMPLETE FOR FLYWAY!!!)
alldat$site %<>% as.character()
alldat$site <- ifelse(
  alldat$site %in% site_lookup$oldsite, 
  site_lookup$mysite[match(alldat$site, site_lookup$oldsite)], 
  alldat$site)

alldat %<>% mutate(latitude=as.numeric(latitude), 
                   longitude=as.numeric(longitude)) %>% 
  filter(metal != "2510.2109999999998")

## remove rows w/ no info in any column
# nas <- is.na(alldat[,-1])
# tormv <- which(rowSums(nas) == ncol(nas))
# alldat <- alldat[-tormv,]

## bespoke fix of some bad lat coords
alldat$latitude <- ifelse(alldat$latitude > 90, alldat$latitude/1000000, alldat$latitude)

## indicate these are color-ring data
alldat$obssource <- rep("cr")

# alldat$bird_id <- ifelse( # combine cr and metal code, use _ 
#   is.na(alldat$combo),    # to signify where one or other is missing
#   alldat$metal, 
#   ifelse(
#     is.na(alldat$metal), 
#     alldat$combo, 
#     paste(alldat$combo, alldat$metal, sep="_")
#   )
# )

# alldat %<>%
#   dplyr::select(bird_id, date, combo, metal, obstype, age, site, latitude, 
#                 longitude, everything())
alldat %<>%
  dplyr::select(date, combo, metal, obstype, age, site, latitude, 
                longitude, everything())

## combine w/ dutch data
dutch <- readRDS("data/analysis/ringing/neth_TPiersma_clean.rds")
alldat <- alldat %>% bind_rows(dutch)

## combine w/ Extremadura/Andalucia data
extr <- readRDS("data/analysis/ringing/spain_JMasero_clean.rds")
alldat <- alldat %>% bind_rows(extr)

## combine w/ French data
fren <- readRDS("data/analysis/ringing/france_FRobin_clean.rds")
alldat %<>% bind_rows(fren)

# combine w/ German - SH data
gsh <- readRDS("data/analysis/ringing/germany_VSalewski_clean.rds")
alldat %<>% bind_rows(gsh)

## combine w/ German - Dummer data
gdu <- readRDS("data/analysis/ringing/germany_JMelter_clean.rds")
alldat %<>% bind_rows(gdu)

## combine w/ Irish - Cork Harbour data
ire <- readRDS("data/analysis/ringing/ireland_BOmahony_clean.rds")
alldat %<>% bind_rows(ire)

##
dutch2 <- readRDS("data/analysis/ringing/neth_MRoodbergen_clean.rds")
alldat %<>% bind_rows(dutch2)


## SAVE ## --------------------------------------------------------------------

# saveRDS(alldat, "data/analysis/ringing/cr_tag_ice_wash_merge.rds")
saveRDS(alldat, "data/analysis/ringing/cr_merge.rds")
