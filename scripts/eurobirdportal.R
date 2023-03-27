## Eurobirdportal data

pacman::p_load(sf, dplyr, magrittr, mapview, stringr, ggplot2)

## centers of 10x10km grid cells

cc <- read.delim("data/citsci/eurobirdportal/cell_centers_10x10/xy_10x10_original.csv", sep=";", header = T)

cc$lon <- as.numeric(scan(text=cc$x, dec=",", sep="."))
cc$lat <- as.numeric(scan(text=cc$y, dec=",", sep="."))

cc %<>% dplyr::select(-x,-y)

cc %>% filter(!is.na(lat)) %>% 
  st_as_sf(coords = c("lon", "lat"), 
                crs = 4326, agr = "constant") %>% 
  slice(165000:180000) %>% 
  mapview()

## filter to cells w/in study area
ccbbox <- cc %>% filter((lat > 0 & lat < 71.5) & (lon > -30 & lon < 45))

ccbbox <- cc


## checklist data 

## fixed and complete checklists
comp_chck <- read.table(
  "data/citsci/eurobirdportal/2023_03_13_Limosa/limlim_lists/limlim_lists.csv",
  sep = ";", header = T)

## casual checklists
cas_chck <- read.table(
  "data/citsci/eurobirdportal/2023_03_13_Limosa/limlim_casual/limlim_casual.csv",
  sep = ";", header = T)

## combine all checklists

chck <- bind_rows(comp_chck, cas_chck)

## only checklists w/in study area
chckbbox <- chck %>% 
  rename(CELLCODE10 = code10x10) %>% 
  filter(CELLCODE10 %in% ccbbox$CELLCODE10)

## summarise n checklists per cell (i.e., events w/ detection)
cell_summ <- chckbbox %>% group_by(CELLCODE10) %>% 
  summarise(
    n_events   = sum(events),
    n_species_records = sum(records_of_species),
    n_observers = sum(observers),
    n_dates    = n_distinct(event_date),
    early_date = min(as.Date(event_date)),
    late_date  = max(as.Date(event_date))
  )

## add cell-center coords to locations w/ godwit records

cell_summ2 <- cell_summ %>% left_join(ccbbox)

cell_summ2_sf <- cell_summ2 %>% 
  st_as_sf(coords = c("lon", "lat"), 
         crs = 4326, agr = "constant") 

cell_summ2_sf %>% 
  mapview(zcol="n_events")


### map it

library(rworldmap)
# get world map
wmap <- getMap(resolution="high")

## Lambert EA: EPSG:3035 ## Lambert CConic: EPSG:3034
wmap_prj <- st_as_sf(wmap) %>% st_transform(crs = "EPSG:3035")

# grid_prj <- globgrid %>% st_transform(crs = "EPSG:3035")
bbox_prj <- st_bbox(
  c(xmin = -16, xmax = 21,
    ymin = 8, ymax = 65.2), crs = 4326) %>% 
  st_as_sfc() %>% st_transform("EPSG:3035") %>% st_bbox()


map <- ggplot() +
  geom_sf(data = wmap_prj, fill = "grey70", color = "white", size=.1) +
  geom_sf(data = arrange(cell_summ2_sf, n_events),
          aes(color = n_events), size=.5) +
  # geom_sf(data = wmap_prj, fill = NA, color = "white", size=0.2) +
  coord_sf(xlim =
             c(bbox_prj[1], bbox_prj[3]), ylim = c(bbox_prj[2], bbox_prj[4]),
           expand = T) +
  viridis::scale_color_viridis(option = "D") +
  theme(
    plot.background = element_rect(fill = "white"),
    panel.background = element_rect(fill = "white"),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    panel.border = element_rect(colour = "black", fill=NA, size=.5),
    legend.key=element_blank()
  ) +
  labs(fill='Checklists') 

ggsave(paste0("figures/10x10_eurobirdportal_nevents.png"), plot=map, width=5, height = 6)
