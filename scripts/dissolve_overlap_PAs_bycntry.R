## So far unsuccessful attempt at dissolving overlapping protected areas 
# (while keeping non-overlaps as separate polygons, and keeping some attr info)

pacman::p_load(dplyr, sf, mapview, wdpar, stringr)

## Country - EEZ union
cntryeez <- st_read("C:/Users/Martim Bill/Documents/political_connectivity/spatial_data/shapefiles_EEZ_countries/union_countries_EEZs/EEZ_Land_v3_202030.shp")


#### Making use of the wdpar package

## country names
x3 <- readRDS("data/analysis/protectedness/wdpa_EAF_ISO.rds")

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
  
  ## clip PAs by country-EEZ borders
  
  ceez <- filter(cntryeez, ISO_TER1 == x3[i]) %>%
    st_transform(crs=st_crs(diss))
  
  diss_int <- st_intersection(diss, ceez) %>% dplyr::select(colnames(diss))
  
  ## save
  saveRDS(
    diss_int, 
    paste0("data/geodata/wdpa/pas_by_country_dissolved/", x3[i], ".rds")
    )
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

