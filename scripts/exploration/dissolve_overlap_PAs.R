## So far unsuccessful attempt at dissolving overlapping protected areas 
# (while keeping non-overlaps as separate polygons, and keeping some attr info)

pacman::p_load(dplyr, sf, mapview)

one <- sf::st_read("data/geodata/wdpa/pas_by_country/PARENT_ISO_ALB.gpkg")
# mapview(one)

## union

st_combine(one) %>%
  st_make_valid()

st_union(one, by_feature = TRUE)

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



