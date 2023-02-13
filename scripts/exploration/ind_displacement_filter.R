### Filter individual data based on displacement threshold 

pacman::p_load(
  dplyr, igraph, stringr, tictoc, tidygraph, sfnetworks, ggplot2, 
  sf, mapview, magrittr, lubridate, netrankr, dggridR, data.table)


## Run through each data type -------------------------------------------------
datatype <- "metal"
# datatype <- "color"
# datatype <- "trax"

if(datatype == "color"){
  ## color ringed bird captures and resightings overlaid on polygon layer
  alldat <- readRDS("data/analysis/ringing/cring_merge_no7dayreobs_ibas.rds")
  alldat %<>% rename(timestamp = date)
} else if (datatype == "metal"){
  ## metal ring captures, recaptures, recoveries
  alldat <- readRDS("data/analysis/ringing/euring_metal_ibas.rds")
  alldat %<>% rename(timestamp = date)
} else if (datatype == "trax"){
  ## tracking locations 
  # alldat <- readRDS("data/analysis/tracking/PTT_GPS_mconn_12h_ibas.rds")
  alldat <- readRDS("data/analysis/tracking/PTT_GPS_mconn_12h_no0_ibas.rds")
  # alldat %<>% rename(timestamp = date)
}

# fxn for splitting string into columns
source("C:/Users/Martim Bill/Documents/R/source_scripts/str2col.R")
# fxn for calc. displacement distance from first tracking locations
source("C:/Users/Martim Bill/Documents/R/source_scripts/displ_dist.R")


## choose data subset to run --------------------------------------------------

## which network to create
season <- "all"
# season <- "spring"
# season <- "fall"

### Separate data for fall and spring migration
## Spring: January 1 - June 30th, Fall: June 24th - January 31st

if(season == "all"){ ## year-round
  netdat <- alldat ## all individuals
} else if(season == "spring"){
  netdat <- subset(alldat, month(alldat$timestamp) %in% c(1:6))
} else if(season == "fall"){
  doy <- lubridate::yday(alldat$timestamp) # June 23/24: 175
  netdat <- subset(
    alldat, 
    month(alldat$timestamp) %in% c(1, 7:12) | doy %in% c(175:181)
  )
}

## Tracking data: rmv ids w/ only local displacement in season (i.e. no mig.)
netdat <- rename(netdat, id = bird_id)

## loop thru IDs adding column of diplacement distance from first location 
tictoc::tic()
netdat <- displ_dist(netdat)
tictoc::toc()

## Summarise max/avg distance from first points -----------------------------

displ_id <- netdat %>% group_by(id) %>% 
  summarise(
    mn_displ = mean(disp_km),
    sd_displ = sd(disp_km),
    md_displ = median(disp_km),
    mx_displ = max(disp_km)
  )

write.csv(displ_id, paste0("data/analysis/summaries/", datatype,"_id_displ_summ.csv"))

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
  geom_vline(aes(xintercept=200), color = "red")

ggsave(paste0("figures/", datatype,"_id_displ_hist.png"), width = 5, height = 4)


## remove individuals w/ only local displacement (i.e. no migration) ----------

wmig <- filter(displ_id, mx_displ >= 100)

allids <- n_distinct(netdat$id)
netdat %<>% filter(id %in% wmig$id)
(1 - n_distinct(netdat$id) / allids) * 100

netdat <- rename(netdat, bird_id = id)


### SAVE ###

if(datatype == "color"){
  saveRDS(netdat, paste0("data/analysis/tracking/cring_merge_no7dayreobs_ibas_", season, ".rds"))
  
} else if (datatype == "metal"){
  saveRDS(netdat, paste0("data/analysis/ringing/euring_metal_ibas_", season, ".rds"))
  
} else if(datatype == "trax"){
  saveRDS(netdat, paste0("data/analysis/tracking/PTT_GPS_mconn_12h_no0_ibas_", season, ".rds"))
}
