## So far unsuccessful attempt at dissolving overlapping protected areas 
# (while keeping non-overlaps as separate polygons, and keeping some attr info)

pacman::p_load(dplyr, sf, mapview, wdpar, stringr)

one <- sf::st_read("data/geodata/wdpa/pas_by_country/PARENT_ISO_ALB.gpkg")
# mapview(one)

## union

comb <- st_combine(one) %>%
  st_make_valid()
mapview(comb)

pa_union <- st_union(one, by_feature = F)
mapview(pa_union)

pa_union2 <- st_union(one, by_feature = TRUE)
mapview(pa_union2)


x <- one %>%
  dplyr::group_by(WDPAID) %>% 
  # dplyr::summarise(across(geometry, ~ sf::st_combine(.)), .groups = "keep") %>%
  # dplyr::summarise(across(geometry, ~ sf::st_union(.)), .groups = "drop")
  dplyr::summarise()

parts <- one %>% st_union() %>% 
  st_cast("POLYGON") %>% st_as_sf()
mapview(parts)

parts <- parts %>% mutate(
  part_id = do.call(
    rbind, 
    lapply(st_intersects(parts, one), function(x) x[1])
  )
)

one$part_id <- 1:nrow(one)

parts <- left_join(parts, st_drop_geometry(one), by = "part_id")


#### Making use of the wdpar package

## country names
folder <- "data/geodata/wdpa/pas_by_country/"

files <- list.files(folder, full.names = F)

x1 <- str_remove(files, pattern = "PARENT_ISO_")
x2 <- str_remove(x1, pattern = ".gpkg")
x3 <- x2[which(!str_detect(x2, pattern = ";"))]

for(i in seq_along(x3)){
  print(x3[i])
  raw <- wdpa_fetch(
    x3[i], wait = FALSE, download_dir = rappdirs::user_data_dir("wdpar")
  )
  
  # clean data
  clean <- wdpa_clean(
    raw,
    exclude_unesco = TRUE, # no biosphere reserves
    retain_status  = c("Designated", "Inscribed", "Established"),
    erase_overlaps = FALSE) 
  
  if(nrow(clean) == 0) {next}
  
  ## remove marine only PAs
  
  clean <- filter(clean, MARINE != "marine")
  
  diss <- wdpa_dissolve(clean)
  
  ## save
  saveRDS(diss, paste0("data/geodata/wdpa/pas_by_country_dissolved/", x3[i], ".rds"))
}


## load and combine all country-dissolved PAs ---------------------------------

folder <- "data/geodata/wdpa/pas_by_country_dissolved/"

files <- list.files(folder, full.names = T)

allpas <- bind_rows(
    lapply(
    seq_along(files), function(x) {
      one <- readRDS(files[x])
    }
  )
)

pa_union <- st_union(allpas, by_feature = F)

