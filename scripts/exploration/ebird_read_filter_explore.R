## Exploring reading in Ebird obs and sampling data, gridding and mapping it

pacman::p_load(auk, mapview, sf, dplyr, terra)
# path to the ebird data file, here a sample included in the package
# get the path to the example data included in the package
# in practice, provide path to ebd, e.g. f_in <- "data/ebd_relFeb-2018.txt
f_in <- "data/ebird/ebd-datafile-SAMPLE_JUL22/ebd_US-AL-101_202204_202204_relApr-2022_SAMPLE/ebd_US-AL-101_202204_202204_relApr-2022.txt"
# output text file

f_out <- "data/ebird/filtered/filtered.txt"

## read in full dataset
ebird_data <- read_ebd(f_in)

## spatialize
# ebird_data %>% 
#   select(longitude, latitude) %>% 
#   st_as_sf(coords = c("longitude", "latitude"), 
#            crs = 4326, agr = "constant") %>% mapview()

# # filter input
# ebird_data_f <- f_in %>%
# # 1. reference file
# auk_ebd() %>%
# # 2. define filters
# auk_species(species = "Black-tailed Godwit") %>%
# # auk_country(country = "Portugal") %>%
# # 3. run filtering
# auk_filter(file = f_out, overwrite = T) %>%
# # 4. read text file into r data frame
# read_ebd()

## spatialize
ebird_data_f %>% 
  select(longitude, latitude) %>% 
  st_as_sf(coords = c("longitude", "latitude"), 
           crs = 4326, agr = "constant") %>% mapview()


## filter by bounding box

boxll <- c(-16.6, 20.48, -16.375 ,20.68)

boxll <- st_bbox(c(xmin = -86.51, xmax = -86.41, 
                   ymin = 32.33, ymax = 32.43),
                 crs = 4326)

boxll %>% st_as_sfc() %>% mapview()

## filter input
ebird_data_f2 <- f_in %>% 
  auk_ebd() %>%
  # 1. spatial bounding box
  auk_bbox(bbox=boxll) %>% 
  auk_filter(file = f_out, overwrite = T) %>% 
  read_ebd()

x <- ebird_data_f2 %>% 
  select(longitude, latitude) %>% 
  st_as_sf(coords = c("longitude", "latitude"), 
           crs = 4326, agr = "constant")
x2 <- boxll %>% st_as_sfc()

mapview(x2) + mapview(x)

### Create grid of checklists per cell ----------------------------------------

## read in full dataset
ebird_data <- read_ebd(f_in)

table(ebird_data$approved)

ebird_data <- ebird_data %>% 
  select("observer_id", "group_identifier", "sampling_event_identifier", 
         "scientific_name", "observation_count", "longitude", "latitude",
         "observation_date")

ebd_summ <- ebird_data %>% 
  group_by(sampling_event_identifier) %>% 
  summarise(
    longitude = first(longitude),
    latitude  = first(latitude)
  )

## Grid the checklists
ebd_vect <- ebd_summ %>% 
  st_as_sf(coords = c("longitude", "latitude"), 
           crs = 4326, agr = "constant") %>% 
  vect()

xtnt <- ext(ebd_vect)

sarea <- rast(xtnt)
# dim(sarea) <- c(5, 5) # grid dimensions

arast <- terra::rasterize(ebd_vect, sarea, fun=sum, background = 0)

x <- raster::raster(arast)
mapview(x)

### Real data -----------------------------------------------------------------
f_in <- "data/ebird/ebd_GB_relJun-2022/ebd_GB_relJun-2022.txt"

ebd <- auk_ebd(f_in)
cols <- c("observer_id", "group_identifier", "sampling_event_identifier", "scientific_name", 
          "observation_count", "longitude", "latitude", "observation_date")

sel <- auk_select(x=ebd, 
                  select = cols,
                  file = f_out, overwrite = T)
xx <- read_ebd(sel)


### Sampling events data (no species records) ---------------------------------

f_in  <- "data/ebird/ebd_sampling_relJun-2022/ebd_sampling_relJun-2022.txt"
# f_in  <- system.file("extdata/zerofill-ex_sampling.txt", package = "auk") # example data
f_out <- "data/ebird/filtered/filter_sampling_relJun-2022.txt"

sam <- auk::auk_sampling(f_in)

## filter by bounding box
boxll <- st_bbox(
  c(xmin = -25.3, xmax = 13.58, # roughly portugal
    ymin = 25.1, ymax = 66.9),
  crs = 4326
)

# boxll %>% st_as_sfc() %>% mapview()

## Remove unneeded columns filter by bounding box
cols <- c("sampling event identifier", "longitude", "latitude", "observer_id",
          "observation date", "number observers", "group identifier")

tictoc::tic()
sel <- f_in %>% 
  auk_sampling() %>%
  # 1. spatial bounding box
  auk_bbox(bbox=boxll) %>% 
  auk_filter(file = f_out, 
             keep = cols,
             # drop = cols,
             overwrite = T)
tictoc::toc()

tictoc::tic()
sam <- auk::read_sampling(sel)
tictoc::toc()


## Grid the checklists
ebd_vect <- sam %>% 
  st_as_sf(coords = c("longitude", "latitude"), 
           crs = 4326, agr = "constant") %>% 
  vect()

xtnt <- ext(ebd_vect)

sarea <- rast(xtnt)
dim(sarea) <- c(100, 50) # grid dimensions

arast <- terra::rasterize(ebd_vect, sarea, fun=sum, background = NA)
plot(arast)

x <- raster::raster(arast)
mapview(x)
