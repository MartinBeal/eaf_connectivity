## processing WDPA data 

pacman::p_load(wdpar)

# download protected area data for Portugal
# (excluding areas represented as point localities)
port_raw <- wdpa_fetch(
  "Portugal", wait = TRUE, download_dir = rappdirs::user_data_dir("wdpar")
)

# clean data
port_clean <- wdpa_clean(
  port_raw,
  exclude_unesco = TRUE, # no biosphere reserves
  retain_status  = c("Designated", "Inscribed", "Established"),
  erase_overlaps = FALSE) 

## remove marine only PAs

port_clean <- filter(port_clean, MARINE != "marine")

port_dissolve <- wdpa_dissolve(port_clean)

port_split <- port_dissolve %>% 
  st_cast("POLYGON") %>% st_as_sf()

## just sf (no wdpdar)
xx <- port_clean %>% st_union() %>% 
  st_cast("POLYGON") %>% st_as_sf()

par(mfrow = c(1, 2))
plot(sf::st_geometry(port_split),
     main = "original", col = "transparent")
plot(sf::st_geometry(xx),
     main = "dissolved", col = "transparent")

port_reid <- port_split %>% mutate(
  part_id = do.call(
    rbind, 
    lapply(st_intersects(port_split, port_clean), function(x) x[1])
  )
)

port_clean$part_id <- 1:nrow(port_clean)

port_reid <- left_join(port_reid, st_drop_geometry(port_clean), by = "part_id")
