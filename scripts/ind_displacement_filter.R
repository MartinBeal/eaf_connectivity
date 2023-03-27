# -----------------------------------------------------------------------------
### Combine all three data type 
## Filter individual data based on displacement threshold -------------------
# -----------------------------------------------------------------------------

pacman::p_load(
  dplyr, igraph, stringr, tictoc, tidygraph, sfnetworks, ggplot2, 
  sf, mapview, magrittr, lubridate, netrankr, dggridR, data.table)


## ring relocations overlaid on polygon layer
coldat <- readRDS("data/analysis/ringing/cring_merge_no7dayreobs.rds")
coldat %<>% rename(timestamp = date) %>% mutate(
  datatype = "color",
  device   = NA,
  study_site = NA
  ) %>% dplyr::select(
    bird_id, timestamp, latitude, longitude, datatype, metal, scheme_country, 
    study_site, device, obstype
  )

## metal ring locations overlaid on polygon layer
metdat <- readRDS("data/analysis/ringing/euring_filtbycr.rds") 
metdat %<>% rename(timestamp = date) %>% mutate(
  datatype="metal",
  device   = NA,
  study_site = NA) %>% dplyr::select(
    bird_id, timestamp, latitude, longitude, datatype, metal, scheme_country, 
    study_site, device, obstype
  )

## tracking locations overlaid on polygon layer (OLD)
# trxdat <- readRDS("data/analysis/tracking/PTT_GPS_mconn_12h_no0_stpovrs_ibas.rds")
## tracking locations speed filtered, resampled
trxdat <- readRDS("data/analysis/tracking/PTT_GPS_mconn_12h_no0.rds")

trxdat %<>% rename(
  metal = animal_ring_id,
  bird_id = id) %>% mutate(
  datatype  = "trax",
  scheme_country = NA,
  obstype = NA,
  ) %>% dplyr::select(
    bird_id, timestamp, latitude, longitude, datatype, metal, scheme_country, 
    study_site, device, obstype
  )

## combine datatype sets ------------------------------------------------------
alldat <- bind_rows(
  coldat, metdat, trxdat
)

# fxn for splitting string into columns
source("C:/Users/Martim Bill/Documents/R/source_scripts/str2col.R")
# fxn for calc. displacement distance from first tracking locations
source("C:/Users/Martim Bill/Documents/R/source_scripts/displ_dist.R")


## choose data subset to run --------------------------------------------------

## Tracking data: rmv ids w/ only local displacement in season (i.e. no mig.)
alldat <- rename(alldat, id = bird_id)

## loop thru IDs adding column of diplacement distance from first location 
tictoc::tic()
alldat <- displ_dist(alldat)
tictoc::toc()


## Summarise max/avg distance from first points -----------------------------

displ_id <- alldat %>% group_by(id) %>% 
  summarise(
    mn_displ = mean(disp_km),
    sd_displ = sd(disp_km),
    md_displ = median(disp_km),
    mx_displ = max(disp_km)
  )

displ_id <- alldat %>% group_by(id, datatype) %>% 
  summarise(
    mn_displ = mean(disp_km),
    sd_displ = sd(disp_km),
    md_displ = median(disp_km),
    mx_displ = max(disp_km)
  )

write.csv(displ_id, paste0("data/analysis/summaries/alldtypes_id_displ_summ.csv"))

## histogram
ggplot(displ_id, mapping=aes(mx_displ)) + 
  xlim(-100, 6500) +
  # ylim(0,30) + # tracking 
  # ylim(0,75) + # metal
  # ylim(0,2050) + # color
  geom_histogram(bins=100) + 
  theme_bw() +
  theme(
    axis.text = element_text(size=11)
  ) + 
  geom_vline(aes(xintercept=100), color = "blue") + 
  # geom_vline(aes(xintercept=200), color = "red") +
  facet_wrap(~datatype, scales="free_y", nrow=3, ncol=1)

ggsave(paste0("figures/alldtypes_id_displ_hist.png"), height = 10, width = 4)


## remove individuals w/ only local displacement (i.e. no migration) ----------

wmig <- filter(displ_id, mx_displ >= 100)

allids <- n_distinct(alldat$id)
alldat2 <- alldat %>% filter(id %in% wmig$id)
(1 - n_distinct(alldat2$id) / allids) * 100 ## % of all IDs removed by filter

alldat2 <- rename(alldat2, bird_id = id)


### SAVE ### ------------------------------------------------------------------

saveRDS(alldat2, "data/analysis/combined/alldatatypes_100km_all.rds")
