### Load extremadura capture, recapture and resighting database - clean, combine

## load packages
pacman::p_load(dplyr, magrittr, lubridate, stringr)

extr_s <- readxl::read_xlsx('data/color/JMasero_spain/Recuperaciones.xlsx')
extr_c <- readxl::read_xlsx('data/color/JMasero_spain/Recuperaciones.xlsx', sheet = 2)

## fix commas in lat lons

extr_c %<>%
  mutate(
    Latitude  = as.numeric(str_replace(Latitude, pattern = ',', '.')),
    Longitude = as.numeric(str_replace(Longitude, pattern = ',', '.'))
    )

## observation type
extr_c$obstype <- "N"
extr_s$obstype <- "R/S" # recaps and resights not distinguished 

## remove birds ringed under Dutch Scheme (duplicated in Jos/Theunis dataset)
dutchids <- unique(subset(extr_c, Colour_ring_code == "Dutch Scheme")$Ring)

extr_s <- subset(extr_s, !Ring %in% dutchids)
extr_c <- subset(extr_c, !Ring %in% dutchids)

## remove birds ringed in Portugal (??duplicated in Zés dataset??)
# portids <- unique(subset(extr_c, Capture_site == "Giganta ricefields. Vila Franca de Xira. Lisboa. Portugal")$Ring)
# 
# extr_s <- subset(extr_s, !Ring %in% portids)
# extr_c <- subset(extr_c, !Ring %in% portids)

## get list of capture sites (for later mapping) ------------------------------

capsites <- extr_c %>%
  group_by(Capture_site) %>% 
  dplyr::summarise(
    latitude = first(Latitude),
    longitude = first(Longitude)
  )

# write.csv(capsites, "data/color/JMasero_spain/capsites.csv", row.names = F)

## fill in combos and age from cap data to resight data -----------------------
agecombos <- extr_c %>%
  group_by(Ring) %>% 
  dplyr::summarise(
    age = first(Age),
    combo = first(Colour_ring_code)
  )

extr_s %<>% left_join(agecombos)

## combine sightings/recaps w/ original ringing caps --------------------------
extr_c %<>% 
  dplyr::rename(
    date = Capture_date, site = Capture_site, metal = Ring, age = Age, 
    combo = Colour_ring_code, latitude = Latitude, longitude = Longitude)

extr_s %<>% 
  dplyr::rename(
    date = Resighting_date, site = Site, metal = Ring, latitude = Latitude, 
    longitude = Longitude) %>% 
  mutate(
    latitude  = as.numeric(str_replace(latitude, pattern = ',', '.')),
    longitude = as.numeric(str_replace(longitude, pattern = ',', '.'))) %>% 
  dplyr::select(-Observer, -Observations, -Capture_date, -Species)

## combine
extr_all <- extr_c %>% dplyr::select(
  metal, combo, date, obstype, site, age, latitude, longitude) %>% 
  bind_rows(extr_s) %>% 
  arrange(metal, date)

extr_all <- extr_all %>% 
  mutate(
    # bird_id = paste(combo, metal, sep="_"),
    area = NA, county = NA, region = NA, country = NA,
    scheme = "spain_JMasero",
    obssource = "cr") %>% 
  select(
    date, combo, metal, obstype, age, site, latitude, longitude, area,
    county, region, country, scheme, obssource)

## Convert age codes to juvenile/adult ----------------------------------------

## SAVE ## --------------------------------------------------------------------

saveRDS(extr_all, "data/analysis/ringing/spain_JMasero_clean.rds")
