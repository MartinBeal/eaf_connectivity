## EURING data prep - BTGO

pacman::p_load(dplyr, magrittr, stringr, lubridate)

all <- read.delim("data/EURING/20220607_Alves/20220607_Alves.psv", sep = "|")

spp <- read.delim("data/EURING/ECSpeciesPipeDelimited.psv", sep = "|")
sch <- read.delim("data/EURING/ECSchemesPipeDelimited.psv", sep = "|")

unique(all$sp)

all <- left_join(all, spp[, 1:2], by=c("SpeciesS"="Code"))
all %<>% dplyr::rename(scientific_name = Name)
table(all$scientific_name)


sch %<>% mutate(Code = str_remove(Code, "  "))
all <- left_join(all, sch[, 1:2], by=c("Scheme"="Code"))
all %<>% dplyr::rename(scheme_country = Country)


### Limosa limosa -------------------------------------------------------------
ll <- all[which(str_detect(all$scientific_name, "Limosa limosa")), ]
table(ll$scheme_country)
rm(all)

## add leading zeroes back in where needed
ll$Date <- ifelse(
  nchar(trunc(ll$Date)) == 7,
  paste(0, ll$Date, sep = ""),
  ll$Date
  )

ll$Date <- as.Date(as.character(ll$Date), "%d%m%Y")

## Map
# sf::st_as_sf(ll, 
#              coords = c("Longitude", "Latitude"), 
#              crs = 4326, agr = "constant") %>% 
#   mapview::mapview()

## Clean up dataframe ---------------------------------------------------------

## Remove birds w/out a recognizable scheme
ll <- subset(ll, !is.na(scheme_country))

## Remove "ringing only" records (no resighting) 
ll <- subset(ll, Ring %in% names(which(sort(table(ll$Ring)) != 1)))

## only columns helping to understand context of obs
ctxt <- ll %>% dplyr::select(c(3,21,9,10:11,17,27,28,29,62))

## retain relevant columns
# ll %<>% dplyr::select(c(1:3, 4:5, 8:11, 14, 16:17, 21:22, 24, 26:29, 57, 59:62))

# Status/obsmethod: New (ie, ringed for the first time), Recapture, Sighting or X=dead)
# X - dead given if the bird was found dead, very sick and not released, or put in captivity and not released
ll$obstype <- ifelse(
  ll$Cmethod == '-' & ll$Circumstance %in% c(28, 80:89), "S", # Resighting (no capture)
  ifelse( # Found dead, killed or caged and not released
    ll$Cmethod == '-' &
      (ll$Condition %in% c(1:3, 5:6) | 
         ll$Circumstance %in% c(8, 10:16, 19:26)), 
    "X", NA)
  )

# Circumstances including death
ifelse(ll$Cmethod == '-' & ll$Circumstance %in% c(8, 10:16, 19:26), "X", NA) # Dead (either found or killed)

## 

ll %<>% dplyr::rename(latitude = Latitude, longitude = Longitude)

ll$obssource <- rep("euring") # distinguish from cr resighting data from schemes

## remove dots from ring #s
ll$Ring <- str_remove_all(ll$Ring, pattern = fixed("."))

ll %<>% 
  mutate(
    Ring = str_remove_all(Ring, pattern = fixed(".")),
    bird_id = Ring,
    combo   = NA
  ) %>% 
  dplyr::select(bird_id, Ring, combo, Date, everything()) %>% 
  dplyr::rename(metal = Ring, date = Date) %>% 
  as_tibble()

## 

ll %<>% dplyr::select(
  bird_id, metal, combo, date, latitude, longitude, CoordAcc, obstype, 
  obssource, scheme_country
)

## remove records w/ spatial accuracy > 20 km (~900)
ll %<>% filter(CoordAcc < 4)

## keep only records from 1980-present
ll %<>% filter(year(date) >= 1980)

## remove ring re-sighint records (keeping only recaptures and recoveries)
ll_noresight <- filter(ll, obstype %in% "S")

## separate finnish resightigns (color-ring records)
ll_fin <- filter(ll, scheme_country == "Finland")

## map
ll %>% filter(obstype == "S") %>%
    sf::st_as_sf(coords = c("longitude", "latitude"),
                 crs = 4326, agr = "constant") %>%
    mapview()

## SAVE -----------------------------------------------------------------------

saveRDS(ll, "data/analysis/ringing/EURING_BTGO.rds")
saveRDS(ll_noresight, "data/analysis/ringing/EURING_BTGO_noresight.rds")
saveRDS(ll_fin, "data/analysis/ringing/EURING_BTGO_finland_cr.rds")

