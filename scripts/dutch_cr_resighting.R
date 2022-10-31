### Load dutch capture, recapture and resighting database - clean, combine

## load packages
pacman::p_load(dplyr, magrittr, plyr, lubridate)

neth_c <- readxl::read_xlsx('data/color/TPiersma_netherlands/Captures.xlsx')
neth_s <- readxl::read_xlsx('data/color/TPiersma_netherlands/Sightings.xlsx')

## characterize records into observation type
neth_c$obstype <- ifelse(neth_c$`Recapture?` == 0, "N", # new captures
                         ifelse(neth_c$`Recapture?` == 1, "R", # Recapture
                                "U")) # Unknown
neth_s$obstype <- "S"

## remove resightings w/ low uncertainty
neth_s <- subset(neth_s, `Certainty Ringcode` != 3)

## remove resightings w/ incomplete ring readings (missing letters/numbers)
# Birds who've lost the flag-ring - code starts w/ 'A0' - remove
noflag <- unique(
  neth_s[which(grepl( "A0", neth_s$ColourCode, fixed = TRUE)), ]$ColourCode)

neth_s <- subset(neth_s, !ColourCode %in% noflag)
neth_c <- subset(neth_c, !ColourCode %in% noflag)

## remaining observations w/ A are legitimate, A's added to meet coding standard
# "Y1AAAA" e.g. is a bird who escaped w/ only a flag, no other c-rings
aremains <- unique(
  neth_s[which(grepl( "A", neth_s$ColourCode, fixed = TRUE)), ]$ColourCode)

## rmv special code obs
neth_c <- subset(neth_c, !neth_c$ColourCode %in% c("1GFFFF", "1GMMMM","1GUUUU","1LFFFF","1LMMMM","1LUUUU","L0LLLL","M1FFFF","M1MMMM","M1UUUU","U1FFFF","U1MMMM","U1UUUU","X0XXXX"))
neth_s <- subset(neth_s, !neth_s$ColourCode %in% c("1GFFFF", "1GMMMM","1GUUUU","1LFFFF","1LMMMM","1LUUUU","L0LLLL","M1FFFF","M1MMMM","M1UUUU","U1FFFF","U1MMMM","U1UUUU","X0XXXX"))

## which individuals carried a tracking device of any type?
trckdid <- neth_c %>% group_by(ColourCode, Ringnumber) %>% 
  dplyr::mutate(
    Satellitetransmitter = ifelse(is.na(Satellitetransmitter), FALSE, Satellitetransmitter),
    Geolocator           = ifelse(is.na(Geolocator), FALSE, Geolocator),
    Radiotransmitter     = ifelse(is.na(Radiotransmitter), 0, Radiotransmitter)
    ) %>% 
  dplyr::summarise(
    tracker = ifelse(
      any(Satellitetransmitter == TRUE) | 
        any(Geolocator == TRUE) | any(Radiotransmitter == 1),
      TRUE, FALSE)
  )

## combine capture (and recapture) + sightings datasets
neth_c2 <- neth_c %>% 
  dplyr::select(
    FldDate, ColourCode, Ringnumber, obstype,"Age RUG", LocationID, Location, 
    Latitude, Longitude, Region, Country)

comb <- neth_s %>% 
  dplyr::select(
    FldDate, ColourCode, Ringnumber, obstype,"Age RUG", LocationID, Location, 
    Latitude, Longitude, Region, Country) %>% 
  bind_rows(neth_c2) %>% 
  arrange(ColourCode, FldDate)

# comb <- right_join(comb, trckdid) # add tracker info

comb_clean <- comb %>% mutate(
  area = NA, county = NA
) %>% 
  dplyr::rename(
    metal=Ringnumber, combo=ColourCode, date=FldDate, site=Location,
    age=`Age RUG`, latitude=Latitude, longitude=Longitude, region=Region,
    country=Country) %>% 
  bind_rows(neth_c2)

## remove obs w/ 0s in Lat/Lon
comb_clean <- subset(comb_clean, latitude != 0 | longitude != 0)

## remove all ids from birds w/ ringing dates b4 2004 (dates area artificial)
ids <- subset(comb_clean, year(comb_clean$date) < 2004)$combo

## make a universal ID for each bird
comb_clean <- comb_clean[-which(
  is.na(comb_clean$combo) & is.na(comb_clean$metal)
  ), ]

# comb_clean$bird_id <- ifelse( # combine cr and metal code 
#   is.na(comb_clean$combo),     
#   comb_clean$metal, 
#   ifelse(
#     is.na(comb_clean$metal), 
#     comb_clean$combo, 
#     paste(comb_clean$combo, comb_clean$metal, sep="_")
#   )
# )

# comb_clean <- comb_clean %>%
#   select(bird_id, everything())

## scheme name
comb_clean$scheme <- "neth_TPiersma"

## indicate these are color-ring data
comb_clean$obssource <- rep("cr")

## prepare to combine w/ other cring resighting records

comb_clean <- comb_clean %>% mutate(
  metal = as.character(comb_clean$metal),
  age = as.character(comb_clean$age)
  )

# comb_clean %<>% dplyr::select(
#   bird_id, date, combo, metal, obstype, age, site, latitude, longitude, area, 
#   county, region, country, scheme, obssource)
comb_clean %<>% dplyr::select(
  date, combo, metal, obstype, age, site, latitude, longitude, area, 
  county, region, country, scheme, obssource)
## fix the coordinates for two sites in Spain and Senegal

comb_clean$longitude <- ifelse(
  comb_clean$site == "Catarroja, Alfafar, Parque Natural de L'Albufera", 
  -0.367, comb_clean$longitude)
comb_clean$longitude <- ifelse(
  comb_clean$site == "Ndioum, Halrar bridge, Lac Diemgologe", 
  -14.6602, comb_clean$longitude)

## one obs w/ lat/lon flipped
comb_clean <- filter(comb_clean, longitude < 50)

comb_clean$latitude <- ifelse(
  comb_clean$site == "Woudrichem, Merwededijk, Groesplaat", 
  comb_clean$latitude + 10, comb_clean$latitude)

comb_clean$latitude <- ifelse(
  comb_clean$site == "Ootmarsum, Ottershagen", 
  comb_clean$latitude - 2, comb_clean$latitude)

comb_clean %>% filter(site == "Ootmarsum, Ottershagen") %>%
   sf::st_as_sf(coords = c("longitude", "latitude"),
                      crs = 4326, agr = "constant") %>%
   mapview()

## SAVE ## --------------------------------------------------------------------

saveRDS(comb_clean, "data/analysis/ringing/neth_TPiersma_clean.rds")

## calculate number of birds ringed in total and per year