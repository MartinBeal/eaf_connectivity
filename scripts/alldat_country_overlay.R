## Get country that points fall w/in

alldat <- readRDS(
  "data/analysis/combined/alldatatypes_ibas_stpovrs_outsites_all.rds")

cntry <- st_as_sf(rnaturalearthdata::countries50)

cntry <- st_make_valid(cntry)

## small subset 
# alldat <- head(alldat, 1000)

alldat_sf <- alldat %>% st_as_sf(
  coords = c("longitude", "latitude"), 
  crs = 4326)

## overlay points on polygons -------------------------------------------------
ov <- sapply(
  st_intersects(alldat_sf, cntry),
  function(x){
    if(length(x) == 0){x <- 'none'}
    return(x[1])
  })

cntry$rowid <- as.character(1:nrow(cntry))

## combine site info with overlap result
pntscntry <- left_join(data.frame(rowid = ov), st_drop_geometry(cntry))

## combine overlap result back into bird locations w/ site info
alldat <- bind_cols(alldat, pntscntry[,c("name", "iso_a3")])

## name points outside polygons as 'none'
alldat %<>% mutate(
  name = ifelse(is.na(name), "none", name),
  iso_a3 = ifelse(is.na(iso_a3), "none", iso_a3)
)

# alldat_sf <- alldat %>% st_as_sf(
#   coords = c("longitude", "latitude"), 
#   crs = 4326)
# mapview(subset(iba, SitRecID %in% unique(alldat$SitRecID))) + mapview(alldat_sf, zcol="IntName")


## get nearest polygon(s) for points outside ----------------------------------

out <- subset(alldat, name == "none")
# out <- head(out, 1000)

out_sf <- out %>% st_as_sf(
  coords = c("longitude", "latitude"), 
  crs = 4326)

# out_sf %>% mapview(zcol="datatype")

tictoc::tic()
closecntry <- st_nearest_feature(out_sf, cntry)
tictoc::toc()

## combine site info with overlap result
closecntry2 <- left_join(
  data.frame(rowid = as.character(closecntry)), 
  st_drop_geometry(cntry)[, c("rowid", "name", "iso_a3")])

## combine overlap result back into bird locations w/ site info (#3)
out2 <- bind_cols(
  dplyr::select(out, -name, -iso_a3),
  closecntry2[,c("name", "iso_a3")])

out2 %>% st_as_sf(
  coords = c("longitude", "latitude"), 
  crs = 4326) %>% mapview(zcol="name")


# recombine all data and re-arrange ------------------------------------------
## remove outs from all points
all_noout <- subset(alldat, name != "none")

alldat2 <- bind_rows(all_noout, out2) %>% 
  arrange(bird_id, timestamp)

## again name points outside polygons as 'none'
alldat2 %<>% rename(
  country = name
)

## SAVE 
saveRDS(
  alldat2, 
  "data/analysis/combined/alldatatypes_ibas_stpovrs_outsites_all_cntry.rds"
  )
