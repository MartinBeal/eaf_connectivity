### Load Maja R's capture, recapture and resighting database - clean, combine

obs  <- readxl::read_xlsx("data/color/MRoodbergen_netherlands/GodwitObsData.xlsx")
ring <- readxl::read_xlsx("data/color/MRoodbergen_netherlands/GodwitRingData.xlsx")

obs %<>% dplyr::rename(
  date = Datum, combo = Kleurcode, metal = Ringnr, 
  site = Plaatsnaam, latitude = GEOX, longitude = GEOY,
  condition = Conditie
) %>% mutate(age = NA, sex = NA)

ring %<>% dplyr::rename(
  date = Datum, combo = Kleurcode, metal = Ringnr, age = Leeftijd, 
  site = Plaatsnaam, sex = Geslacht, latitude = GEOX, longitude = GEOY
)

unique(obs$metal) %in% unique(ring$metal)

## combine ringing and obs
ring$obstype <- "N"
obs$obstype <- "S"

comb <- ring %>% bind_rows(obs) %>% arrange(metal, combo, date)

## fix a few wayward coordinates

comb$latitude <- ifelse(
  comb$site == "Les Gavines, Racó de l'Olla, PN de L'Albufera, SPANJE",
       39.3, comb$latitude)

## remove individuals that weren't resighted
x <- comb %>% group_by(metal) %>% 
  dplyr::summarise(nobs = n()) %>% filter(nobs == 1)
comb_clean <- filter(comb, !metal %in% x$metal)

## 
comb_clean %<>% filter(!is.na(latitude))
comb_clean %<>% filter(!is.na(longitude))

##
comb_clean$metal <- str_remove(comb_clean$metal, pattern = fixed("."))

## 
comb_clean %<>% mutate(
  obssource = "cr",
  scheme    = "MRoodbergen_netherlands",
  area = NA, county = NA, region = NA, country = NA
)

## 
# comb_clean %>%
#   sf::st_as_sf(coords = c("longitude", "latitude"),
#                crs = 4326, agr = "constant") %>%
#   mapview()


comb_clean %<>% dplyr::select(
  date, combo, metal, obstype, age, site, latitude, longitude, area, 
  county, region, country, scheme, obssource)


## SAVE ## --------------------------------------------------------------------

saveRDS(comb_clean, "data/analysis/ringing/neth_MRoodbergen_clean.rds")
