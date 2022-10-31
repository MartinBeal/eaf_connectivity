### Load irech capture, recapture and resighting database - clean, combine

## load packages
pacman::p_load(dplyr, magrittr, lubridate, stringr)

ire <- readxl::read_xlsx("data/color/BOmahoney_ireland/MB copy database.xlsx", sheet = 2)

ire <- ire[-822,]

ire <- subset(ire, !is.na(ire$`Co-ordinants`))

ire$`Co-ordinants` <- gsub("\"", "", ire$`Co-ordinants`)

ire$latitude <- sapply(str_split(ire$`Co-ordinants`, "N"), function(x) x[[1]])
ire$longitude <- sapply(str_split(ire$`Co-ordinants`, "N"), function(x) x[[2]])

## negative longitude to add in later
neg_lon <- which(str_detect(ire$longitude, "W"))

ire$latitude <- str_trim(
  str_replace_all(ire$latitude, c("W"="", "E"="", "N"=""), ""))
ire$longitude <- str_trim(
  str_replace_all(ire$longitude, c("W"="", "E"="", "N"=""), ""))

hms <- which(substr(ire$`Co-ordinants`, 3, 3) == "°")

lat <- ire$latitude[hms]
lon <- ire$longitude[hms]

lat  <- gsub('([[:punct:]])|\\s+',' ', lat)
lon  <- gsub('([[:punct:]])|\\s+',' ', lon)


# lat <- gsub('([[:punct:]])|\\s+',' ', ire$latitude[hms])
# lon <- gsub('([[:punct:]])|\\s+',' ', ire$longitude[hms])
# 
# lat <- trimws(str_replace_all(lat, c("W"="", "E"="", "N"=""), ""))
# lon <- trimws(str_replace_all(lon, c("W"="", "E"="", "N"=""), ""))


library(measurements)

decdeg_lat <- sapply(seq_along(lat), function(x) {
  one <- lat[x]
  out <- conv_unit(one, from = "deg_min_sec", to = "dec_deg")
  return(out)
} )

decdeg_lon <- sapply(seq_along(lon), function(x) {
  one <- lon[x]
  out <- conv_unit(one, from = "deg_min_sec", to = "dec_deg")
  return(out)
} )

ire$latitude[hms]  <- decdeg_lat
ire$longitude[hms] <- decdeg_lon

ire$latitude  <- gsub("°", "", ire$latitude)
ire$longitude <- gsub("°", "", ire$longitude)

ire$latitude  <-  as.numeric(gsub(" ", ".", ire$latitude))
ire$longitude <-  as.numeric(gsub(" ", ".", ire$longitude))

## add negatives to longitude dec deg
ire$longitude[neg_lon] <- ire$longitude[neg_lon] * -1

ire <- subset(ire, !is.na(latitude))
ire <- subset(ire, !is.na(longitude))

# sf::st_as_sf(ire,
#              coords = c("longitude", "latitude"),
#              crs = 4326, agr = "constant") %>% mapview()


## translate column names 
ire %<>% dplyr::rename(
  date = Date, combo = `Colour-rings`, metal = "Ring number", 
  obstype = `Record type`, site = Location, county = County, 
  country = Country, age = "Age at ringing"
)

##
ire$region <- NA
ire$area <- NA

ire %<>% filter(latitude != 0)

## date to datetime
ire$date <- openxlsx::convertToDate(ire$date)
ire %<>% filter(!is.na(date))

## scheme name
ire$scheme <- "irel_BOmahony"

## indicate these are color-ring data
ire$obssource <- rep("cr")

## observation type 
# ire$obstype <- ifelse(ire$obstype == "s", "S", 
#                        ifelse(ire$obstype == "F", "S", ire$obstype))

## birds missing metal rings
ire$metal <- ifelse(ire$metal == "", NA, ire$metal)

## birds only ringed, never re-seen
x <- ire %>% group_by(metal) %>% dplyr::summarise(n = n()) %>% filter(n==1)
ire <- filter(ire, !metal %in% x$metal)

## get important columns for combining
ire_clean <- ire %>% dplyr::select(
  date, combo, metal, obstype, age, site, latitude, longitude, area,
  county, region, country, scheme, obssource)

## map it
sf::st_as_sf(ire_clean,
             coords = c("longitude", "latitude"),
             crs = 4326, agr = "constant")  %>%
  mapview::mapview()

## SAVE -----------------------------------------------------------------------

saveRDS(ire_clean, "data/analysis/ringing/ireland_BOmahony_clean.rds")
