## International waterbird census data


iwc <- read.csv("data/IWC_counts/DR_MBeal_080223.csv")

## only godwit data

iwc %<>% filter(reporting_name == "Limosa limosa")

iwc <- iwc %>% 
  mutate(
    latitude = y,
    longitude = x
  ) %>% 
  st_as_sf(coords = c("x", "y"), 
         crs = 4326, agr = "constant")

## Remove bad count (wrong species) at Banc D'Arguin
iwc <- filter(iwc, site != "MR00030" & count != 15268)

## Remove locations w/ no 

x <- iwc %>% arrange(desc(count)) %>% slice(1:1000) %>% mapview(zcol="count")

mapshot(x, url="C:/Users/Martim Bill/Desktop/XX.html")

## summarize across years
iwc_summ <- iwc %>% 
  mutate(visit_date = as.Date(visit_date)) %>% 
  group_by(iwccountry, site, sitename) %>%
  arrange(visit_date) %>% 
  summarise(
    n_yrs      = n_distinct(lubridate::year(visit_date)),
    last_year  = lubridate::year(last(visit_date)),
    last_count = last(count),
    mean_count = mean(count),
    sd         = sd(count),
    max_count = max(count),
    min_count = min(count),
    latitude  = median(latitude),
    longitude = median(longitude)
  ) %>% st_as_sf(coords = c("longitude", "latitude"), 
               crs = 4326, agr = "constant")


mapview(iwc_summ, zcol="mean_count")
iwc_summ100 <- iwc_summ %>% ungroup() %>% arrange(desc(mean_count)) %>% slice(1:100) 

iwc_summ100 %>% mapview(zcol="mean_count")


## Overlay on nodes/sites found from tracking/color/metal ---------------------

allsites <- readRDS(
  "data/analysis/site_nodes/allsite_polygons_alldatatypes_all.rds")

allsites$rowid <- as.character(seq_len(nrow(allsites)))

## Load site summary both IBAs and outsite centroids) ------------------------
site_summ <- readRDS(
  paste0("data/analysis/site_nodes/alldatatypes_allsites_cent_all.rds")
)

allsites %<>% left_join(st_drop_geometry(site_summ))

## overlay count points on site polys, and get T/F if overlaps
iwc_summ$ovr_mct_site <- apply(st_intersects(iwc_summ, allsites), 1, any)


## overlay points on site polygons --------------------------------------------
ov <- sapply(
  st_intersects(iwc_summ, allsites),
  function(x){
    if(length(x) == 0){x <- 'none'}
    return(x[1])
  })

allsites$rowid <- as.character(1:nrow(allsites))

## combine site info with overlap result
pntssites <- left_join(data.frame(rowid = ov), st_drop_geometry(allsites))

## combine overlap result back into bird locations w/ site info
iwc_summ <- bind_cols(iwc_summ, pntssites[,c("SitRecID", "site_poly")])

## name points outside polygons as 'none'
iwc_summ %<>% mutate(
  SitRecID = ifelse(is.na(SitRecID), "none", SitRecID),
  site_poly = ifelse(is.na(site_poly), "none", site_poly)
)

## get nearest polygon(s) for points outside ----------------------------------

out <- subset(iwc_summ, SitRecID == "none")

out_sf <- out %>% st_as_sf(
  coords = c("longitude", "latitude"), 
  crs = 4326)

tictoc::tic()
siteinrange <- st_is_within_distance(out_sf, allsites, 10000) ## m from near. poly
tictoc::toc()

## classify result (none, iba row id, multiple)
whichpoly <- sapply(
  siteinrange,
  function(x){
    if(length(x) == 0)
    {x <- 'none'} else if (length(x) > 1)
    {x <- 'multiple'
    }
    return(x[1])
  })

## get rowid of ibas for points falling in just one iba (keep 'nones')
singles_rid <- as.character(whichpoly[which(whichpoly != "multiple")])
singles <- out_sf[which(whichpoly != "multiple"), ] # filter from bird data

## get iba metadata for ibas that points fall within
single_site <- left_join(
  data.frame(rowid = singles_rid), 
  st_drop_geometry(allsites)[, c("rowid", "SitRecID", "site_poly")]) %>% 
  mutate(
    SitRecID = as.character(SitRecID))

## combine overlap result back into bird locations w/ site info (#2)
singles <- bind_cols(
  subset(singles, select=-c(SitRecID, site_poly)),
  single_site[,c("SitRecID", "site_poly")])


## multiple polygons ---------------------------------------------------------
multis   <- out[which(whichpoly == "multiple"), ]
multi_sf <- out_sf[which(whichpoly == "multiple"), ]

## get nearest polygon for those obs w/ multiple within search distance
tictoc::tic()
neariba <- st_nearest_feature(multi_sf, allsites)
tictoc::toc()

## combine site info with overlap result
multis_ibas <- left_join(
  data.frame(rowid = as.character(neariba)), 
  st_drop_geometry(allsites)[, c("rowid", "SitRecID", "site_poly")]) %>% 
  mutate(
    SitRecID = as.character(SitRecID))

## combine overlap result back into bird locations w/ site info (#3)
multis <- bind_cols(
  subset(multis, select=-c(SitRecID, site_poly)),
  multis_ibas[,c("SitRecID", "site_poly")])


## recombine all data and re-arrange ------------------------------------------
## remove outs from all points
all_noout <- subset(iwc_summ, SitRecID != "none")

iwc_summ2 <- bind_rows(all_noout, singles, multis) %>% 
  arrange(site, last_year)

## again name points outside polygons as 'none'
iwc_summ2 %<>% mutate(
  SitRecID = ifelse(is.na(SitRecID), "none", SitRecID),
  site_poly = ifelse(is.na(site_poly), "none", site_poly)
)

### do some sense-checks
## check point falling in/around a site
mapview(subset(allsites, site_poly == "Tejo estuary")) + 
  (iwc_summ2 %>% filter(site_poly == "Tejo estuary") %>% 
     st_as_sf(
       coords = c("longitude", "latitude"),
       crs = 4326) %>% mapview())
# ## points not considered in/around a site 
# mapview(subset(iba, IntName == "Tejo estuary")) + (alldat2 %>% filter(IntName == "none") %>% st_as_sf(
#   coords = c("longitude", "latitude"), 
#   crs = 4326) %>% mapview())


## remove sites w/ 0 counts
iwc_summ2 %<>% filter(min_count > 1)

## SAVE -----------------------------------------------------------------------

saveRDS(
  iwc_summ2,
  "data/analysis/site_counts/IWC_Llimosa_MCTsites_overlay.rds"
)


## maps ------------------------------------------------------------------------

## just count sites and mean totals ------------------------------------------

# get world map
wmap <- rworldmap::getMap(resolution="high")

## Lambert EA: EPSG:3035 ## Lambert CConic: EPSG:3034
# rams_sf_prj <- st_as_sf(rams) %>% st_transform(crs = "EPSG:3035")
# iba_sf_prj <- st_as_sf(iba) %>% st_transform(crs = "EPSG:3035")
wmap_prj <- st_as_sf(wmap) %>% st_transform(crs = "EPSG:3035")

## Map of all sites -----------------------------------------------------------
bbox_prj <- st_bbox( ## EAF bbox
  c(xmin = -16, xmax = 21,
    ymin = 8, ymax = 65.2), crs = 4326) %>% 
  st_as_sfc() %>% st_transform("EPSG:3035") %>% st_bbox()

## save theme to re-use
maptheme <- 
  theme(
    # plot.background = element_rect(fill = "white"),
    panel.background = element_rect(fill = "white"),
    # axis.text = element_blank(),
    # axis.ticks = element_blank(),
    # panel.grid.major = element_line(colour="grey60", size=.2),
    panel.grid = element_blank(),
    panel.border = element_rect(colour = "black", fill=NA, size=.5),
    legend.key=element_blank()
  )

### 
map <- ggplot() +
  geom_sf(data = wmap_prj, fill = "grey70", color = "white", size=0.2) +
  geom_sf(data = iwc_summ2,  
          aes(), color = "black", size=0.5) +
  geom_sf(data = arrange(iwc_summ100, mean_count),  
          aes(color = mean_count), size = 2) +
  coord_sf(xlim = c(bbox_prj[1], bbox_prj[3]), 
           ylim = c(bbox_prj[2], bbox_prj[4]), 
           expand = T) +
  scale_size(range = c(0.01, 2))  +
  scale_color_distiller(palette = "OrRd") +
  maptheme

map

## Save ##
ggsave("figures/IWC_top100_count.png", 
       plot=map, width=5, height = 6)


## Show count sites overlapping with sites IDd from M/C/T data ----------------

iwc_summ2$inout <- ifelse(iwc_summ2$SitRecID == "none", "Out", "In")

map2 <- ggplot() +
  geom_sf(data = wmap_prj, fill = "grey70", color = "white", size=0.2) +
  # geom_sf(data = site_summ,  
  #         aes(), color = "black", size=.5) +
  geom_sf(data = iwc_summ2,  
          aes(color = inout), size=1) +
  coord_sf(xlim = c(bbox_prj[1], bbox_prj[3]), 
           ylim = c(bbox_prj[2], bbox_prj[4]), 
           expand = T) +
  scale_size(range = c(0.01, 2)) +
  maptheme

map2

## Save ##
ggsave("figures/IWC_inout_MCTsites.png", 
       plot=map2, width=5, height = 6)
ggsave("figures/IWC_inout_MCTsites_allsites.png", 
       plot=map2, width=5, height = 6)
