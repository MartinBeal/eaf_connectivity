## Calculate % coverage of BTGO sites by protected areas ----------------------
# -----------------------------------------------------------------------------

pacman::p_load(sf, stringr, dplyr)

source("C:/Users/Martim Bill/Documents/R/source_scripts/str2col.R")


outpolys <- readRDS(
        "data/analysis/site_nodes/outsite_polygons_alldatatypes_all.rds")

## polygon dataset (IBAs) -----------------------------------------------------
iba <- st_read("data/geodata/ibas/EAF_btgo_IBA/eaf_btgo_iba.shp")

iba %<>% st_make_valid() # deal with polygon invalidity

## SAVE 
site_summ <- readRDS(
  "data/analysis/site_nodes/alldatatypes_allsites_cent_all.rds")

iba_visit <- iba %>% filter(SitRecID %in% site_summ$SitRecID)

iba_red <- dplyr::select(
  iba_visit, SitRecID, IntName, Country, Protect, KBACurrent, GISArea) %>% 
  rename(area = GISArea, country = Country)


## join outsite and IBAs
# allsite_poly <- bind_rows(iba_red, outpolys)

## Outsites -  protectedness, area, country -----------------------------------

## area -------------------------------------------------
# NEED TO PROJECT FIRST
outpolys$area <- units::set_units(st_area(outpolys), "km2")

## Countries -------------------------------------------------
## centroids:
out_cent <- ccenters_re %>% 
  mutate(
    cellgrp = paste0("cellgrp_", cellgrp),
    # site_poly = cellgrp, ## make the cell group the main site ID
    SitRecID  = cellgrp ## make the cell group the main site ID
  ) %>% group_by(SitRecID, cellgrp) %>% 
  st_as_sf(coords = c("longitude", "latitude"), 
           crs = 4326, agr = "constant") %>% 
  summarise() %>% 
  dplyr::select(SitRecID, geometry)

cntry <- st_as_sf(rnaturalearthdata::countries50) %>% st_make_valid()

## overlay points on polygons
ov <- sapply(
  st_intersects(out_cent, cntry),
  function(x){
    if(length(x) == 0){x <- 'none'}
    return(x[1])
  })

cntry$rowid <- as.character(1:nrow(cntry))

## combine site info with overlap result
centcntry <- left_join(data.frame(rowid = ov), st_drop_geometry(cntry))

## combine overlap result back into bird locations w/ site info
outpolys <- bind_cols(outpolys, centcntry[,c("sovereignt", "sov_a3")]) %>% 
  rename(country = sovereignt)

## protectedness -------------------------------------------------

## qualitative level
# "some"        "whole"       "little/none" "most"        "unknown"


## ----------------------------------------------------------------------------
## Overlay nodes on PAs polygons per country 
## ----------------------------------------------------------------------------

folder <- "data/geodata/wdpa/pas_by_country_dissolved/"

files     <- list.files(folder, full.names = T)
filenames <- list.files(folder, full.names = F)

cntryeez <- st_read("C:/Users/Martim Bill/Documents/political_connectivity/spatial_data/shapefiles_EEZ_countries/union_countries_EEZs/EEZ_Land_v3_202030.shp")

## polygons of all sites (nodes)
allsites <- readRDS(
  "data/analysis/site_nodes/allsite_polygons_alldatatypes_all.rds")

## split countries in to two columns when site shared by two countries
allsites <- allsites %>% bind_cols(
  str2col(
    allsites$ISO3, 
    pattern = "_", 
    cols = 1:2, 
    colnames=c("country1", "country2"))
)

## ----------------------------------------------------------------------------
## overlay site polygons on country PAs ---------------------------------------

alist <- list()

## still need to run: 31, 63, 78

for(i in seq_along(files)){ #c(31, 63, 78)
  cntrypas <- readRDS(files[i])
  onename <- str_remove(filenames[i], pattern = ".rds")
  print(paste(i, onename))
  
  if(i != 20) {
    if(!st_is_valid(cntrypas)){
      cntrypas <- st_make_valid(cntrypas)
    }
  }
  
  cntrynodes <- filter(
    allsites, 
    allsites$country1 == onename | allsites$country2 == onename) %>% 
    st_transform(crs = st_crs(cntrypas)) %>% 
    mutate(
      GISArea = units::set_units(
        st_area(.), "km2")
    )
  
  if(nrow(cntrynodes) < 1) next
  
  if(any(!st_is_valid(cntrynodes))){
    cntrynodes <- st_make_valid(cntrynodes)
  }
  
  # Calculate area and tidy up
  intersect_pct <- st_intersection(cntrynodes, cntrypas) %>% 
    mutate(intersect_area = 
             units::set_units(
               st_area(.), "km2"),
           pct_cover      =  units::drop_units((intersect_area / GISArea) *100)
    ) %>%   # create new column with shape area
    # dplyr::select(NAME, intersect_area) %>%   # only select columns needed to merge
    st_drop_geometry()  # drop geometry as we don't need it
  
  ## add back in sites w/ no PA overlap
  nointer <- st_drop_geometry(
    cntrynodes[which(!cntrynodes$SitRecID %in% intersect_pct$SitRecID), ]
    ) %>% mutate(
      intersect_area = units::set_units(0, "km2"),
      pct_cover      = 0
    )
  
  intersect_all <- bind_rows(intersect_pct, nointer)
  
  saveRDS(
    intersect_all, 
    paste0("data/analysis/protectedness/bycountry/", onename, ".rds")
    )

}


## Read in and combine country-level PA coverage ------------------------------
folder <- "data/analysis/protectedness/bycountry/"

files     <- list.files(folder, full.names = T)

alist <- lapply(
  seq_along(files), function(x){
    one <- readRDS(files[x])
    # if(any(one$SitRecID == "3586")) print(files[x]) # check certain site
    if(nrow(one) > 0) {return(one)}
  }
)

alist <- alist[!sapply(alist,is.null)]

allsites_cvr <- do.call(rbind, alist)
hist(allsites_cvr$pct_cover)

## for transboundary sites, add up coverage

## which sites are duplicated
dups    <- allsites_cvr[which(duplicated(allsites_cvr$SitRecID)), ]
notdups <- allsites_cvr[which(!duplicated(allsites_cvr$SitRecID)), ]

## rmv dups from only one country
gooddups <- filter(dups, country1 != country2)

allsites_cvr2 <- bind_rows(notdups, gooddups)

allsites_red <- allsites_cvr2 %>% 
  group_by(SitRecID) %>% 
  summarise(
    GISArea = first(GISArea),
    country = paste(unique(country), collapse = "_"),
    country1 = first(country1),
    country2 = first(country2),
    intersect_area = sum(intersect_area),
    pct_cover = (intersect_area / GISArea) * 100
  ) %>% 
  mutate(
    pct_cover = units::drop_units(pct_cover)
  )

allsites_red[which(units::drop_units(allsites_red$pct_cover) > 100.01),]


## summarise by country
cntry_cvr <- allsites_cvr %>% 
  group_by(country1) %>% 
  summarise(
    country = first(country),
    n_sites     = n(),
    n_anycover  = sum(ifelse(as.numeric(pct_cover) > 0, T, F)),
    perc_wcover = (n_anycover / n_sites)*100,
    mn_cover = mean(pct_cover),
    sd_cover = sd(pct_cover),
    mn_site_area  = mean(GISArea),
    sd_site_area  = sd(GISArea)
  )

gen_summ <- allsites_cvr %>% 
  summarise(
    n_sites     = n(),
    n_anycover  = sum(ifelse(as.numeric(pct_cover) > 0, T, F)),
    perc_wcover = (n_anycover / n_sites)*100,
    mn_cover = mean(pct_cover),
    sd_cover = sd(pct_cover),
  )

## SAVE

saveRDS(allsites_red, "data/analysis/protectedness/perc_pa_cover_allsites.rds")
saveRDS(cntry_cvr, "data/analysis/protectedness/perc_pa_cover_countries.rds")
saveRDS(gen_summ, "data/analysis/protectedness/perc_pa_cover_gen.rds")

