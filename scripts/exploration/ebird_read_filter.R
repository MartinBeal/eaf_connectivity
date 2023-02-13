## filter and read ebird sampling events data and map 'effort' ---------------

pacman::p_load(auk, mapview, sf, dplyr, terra, rasterVis, rasterVis, ggplot2)

### Sampling events data (no species records) ---------------------------------

f_in  <- "data/ebird/ebd_sampling_relJun-2022/ebd_sampling_relJun-2022.txt"
# f_in  <- system.file("extdata/zerofill-ex_sampling.txt", package = "auk") # example data
f_out <- "data/ebird/filtered/filter_sampling_relJun-2022.txt"

sam <- auk::auk_sampling(f_in)

## filter by bounding box
boxll <- st_bbox(
  c(xmin = -30, xmax = 35, # roughly the EAF
    ymin = 10, ymax = 70),
  crs = 4326
  )

# boxll %>% st_as_sfc() %>% mapview()

## Remove unneeded columns filter by bounding box
cols <- c("sampling event identifier", "longitude", "latitude", "observer_id",
          "observation date", "number observers", "group identifier")

tictoc::tic()
sele <- f_in %>% 
  auk_sampling() %>%
  # 1. spatial bounding box
  auk_bbox(bbox=boxll) %>% 
  auk_filter(file = f_out, 
             keep = cols,
             # drop = cols,
             overwrite = T)
tictoc::toc()

tictoc::tic()
sam <- auk::read_sampling(sele)
tictoc::toc()


## Grid the checklists

## spatialize and project
ebd_vect <- sam %>% 
  st_as_sf(coords = c("longitude", "latitude"), 
           crs = 4326, agr = "constant") %>% 
  st_transform(crs = "EPSG:3035") %>% # lambert EA Europe crs
  vect()

xtnt <- ext(ebd_vect)

sarea <- rast(xtnt, crs = crs(ebd_vect))
dim(sarea) <- c(100, 80) # grid dimensions

arast <- terra::rasterize(
  ebd_vect, sarea, fun=sum, background = NA)
# plot(arast)

# x <- raster::raster(arast)
# mapview(x)

## map it ---------------------------------------------------------------------
# get world map
library(rworldmap)
wmap <- sf::st_as_sf(getMap(resolution="high"))
wmap_prj <- wmap %>% st_transform(crs = "EPSG:3035")

# arast_prj <- terra::project(arast, "EPSG:3035")
# xtnt_prj <- ext(arast_prj)

map <- rasterVis::gplot(arast) +
  geom_sf(data = wmap_prj, color = NA, fill = "grey65", inherit.aes = FALSE) +
  geom_tile(aes(fill = value)) +
  scale_fill_viridis(
    option = "inferno",
    trans  = "sqrt",
    # trans  ="log2",
    breaks = scales::trans_breaks("sqrt", function(x) x ^ 2),
    # breaks = scales::trans_breaks("log2", function(x) 2 ^ x),
    labels = function(x) round(x, 1),
    na.value = NA,
    direction = -1
  ) +
  scale_color_viridis( # give cell border same color as fill
    option = "inferno",
    trans  = "sqrt",
    # trans  ="log2",
    breaks = scales::trans_breaks("sqrt", function(x) x ^ 2),
    # breaks = scales::trans_breaks("log2", function(x) 2 ^ x),
    labels = function(x) round(x, 1),
    na.value = NA,
    direction = -1
    ) +
  geom_sf(data = wmap_prj, color = "grey15", fill = NA, inherit.aes = FALSE) +
  coord_sf(
    crs = "EPSG:3035",
    xlim = c(xtnt_prj[1], xtnt_prj[2]),
    ylim = c(xtnt_prj[3], xtnt_prj[4]),
    expand = F) +
  theme_bw() +
  theme(
    axis.title = element_blank()
  )
map
ggsave("figures/cit_sci/ebird_nchecklistsX3.png", plot = map, width=5.5, height=5.5)
