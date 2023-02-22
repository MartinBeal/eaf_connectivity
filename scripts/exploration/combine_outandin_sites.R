## Create unified layer of IBAs and outsites used by godwits ------------------

## summary datasets of all sites used by godwits ------------------------------

site_summ <- readRDS(
  "data/analysis/site_nodes/alldatatypes_allsites_cent_all.rds"
  )

## outsites (hexcell groups) --------------------------------------------------
outgrid <- readRDS(
  "data/analysis/site_nodes/outsite_polygons_alldatatypes_all.rds"
  )

## polygon dataset (IBAs) -----------------------------------------------------
iba <- st_read("data/geodata/ibas/EAF_btgo_IBA/eaf_btgo_iba.shp")

iba %<>% st_make_valid() # deal with polygon invalidity

usediba <- filter(iba, SitRecID %in% unique(site_summ$SitRecID))

## combine IBAs and outsites --------------------------------------------------

outgrid %<>% 
  mutate(
    cellgrp = paste0("cellgrp_", cellgrp),
    SitRecID  = cellgrp, ## make the cell group the main site ID
    site_poly = "none",
    protect   = NA
  )

allsites <- usediba %>% 
  dplyr::rename(site_poly = IntName, protect = Protect) %>% 
  dplyr::select(SitRecID, site_poly, protect, geometry) %>% 
  mutate(SitRecID = as.character(SitRecID)) %>% 
  bind_rows(dplyr::select(outgrid, SitRecID, site_poly, protect, geometry))

allsites$area <- units::set_units(st_area(allsites), "km^2")

## save
saveRDS(
  allsites,
  "data/analysis/alldatatypes_allsites_poly_all.rds"
)
