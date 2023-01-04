#### Identify and remove individuals w/out migration data
## HAS BEEN CONVERTED INTO A FUNCTION: displ_dist()

pacman::p_load(dplyr, amt, ggplot2)

# trax <- readRDS("data/analysis/tracking/PTT_GPS_mconn_12h.rds")
trax <- readRDS("data/analysis/tracking/PTT_GPS_mconn_12h_no0.rds") # no 0 locs


## Calculate displacement -----------------------------------------------------
proj <- crsuggest::suggest_crs(trax_sf, limit = 1, drop_na = T)


trax_sf <- st_as_sf(trax, coords = c("longitude", "latitude"), 
                    crs = 4326, agr = "constant")

trax_amt <- trax %>% 
  make_track(.x=longitude, .y=latitude, .t=timestamp, 
             id = id, crs = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_d"), 
             all_cols = T)

trax_amt_proj <- transform_coords(trax_amt, crs_to=as.integer(proj$crs_code))

## loop calculation
alist <- list()
for(i in 1:n_distinct(trax_amt_proj$id)){
  print(i)
  one_proj <- subset(trax_amt_proj, id == unique(trax_amt_proj$id)[i])
  one <- subset(trax_amt, id == unique(trax_amt$id)[i])
  one$nsd <- nsd(one_proj)
  alist[[i]] <- one 
}

trax_amt <- do.call(rbind, alist)

trax$disp_km <- sqrt(trax_amt$nsd) / 1000 # km

## Summarise max/avg distance from first points -------------------------------

displ_id <- trax %>% group_by(id) %>% 
  summarise(
    mn_displ = mean(disp_km),
    sd_displ = sd(disp_km),
    md_displ = median(disp_km),
    mx_displ = max(disp_km)
  )

## remove individuals w/ only local displacement (i.e. no migration)

wmig <- filter(displ_id, mx_displ >= 200)

n_distinct(trax$id)
trax %<>% filter(id %in% wmig$id)
n_distinct(trax$id)


trax %>% filter(id == "R\xfcschke") %>% 
  st_as_sf(coords = c("longitude", "latitude"), 
           crs = 4326, agr = "constant") %>% mapview()

## SAVE ## --------------------------------------------------------------------
saveRDS(trax, "data/analysis/tracking/PTT_GPS_mconn_12h_no0_migids.rds")

