## country-level map of ringing effort of L. limosa

## Ringing totals (from EURING, or from national ringing centers)
rtot <- read.csv("data/EURING/national_ringing_totals.csv")
## Euring scheme codes
# sch <- read.delim("data/EURING/ECSchemesPipeDelimited.psv", sep = "|")
sch <- read.csv("data/EURING/EURING_scheme_codes.csv")
sch$Code <- trimws(sch$Code)

rtot <- left_join(rtot, sch[,c(1,2)], by=c("Scheme"="Code"))

## add country names for old scheme codes (e.g., DK for denmark)
rtot$Country <- ifelse(
  str_detect(rtot$Scheme, pattern="DK"), "Denmark", rtot$Country)

## Quickly country polygons 
library(rworldmap)
# get world map
wmap <- getMap(resolution="high")

## add up ringing totals across years
rtot_all <- rtot %>% 
  group_by(Country) %>% 
  summarise(
    total.ringed = sum(na.omit(Total.ringed)),
    yr_range = paste(min(Year), max(Year), sep="-")
    )

rtot_all$Country <- ifelse(
  rtot_all$Country == "Czech Republic", "Czech Rep.", rtot_all$Country)

## duplicate row from UK & Ireland so that they have same # (schemes combined)
rtot_all$Country <- ifelse(
  rtot_all$Country == "UK & Ireland", "United Kingdom", rtot_all$Country)
x <- subset(rtot_all, Country == "United Kingdom")
x$Country <- x$Country <- "Ireland"
rtot_all <- rbind(rtot_all, x)

tots_sf <- wmap %>% st_as_sf() %>% 
  dplyr::select(NAME) %>% 
  left_join(rtot_all, by = c("NAME"="Country"))
  
mapview(tots_sf, zcol="total.ringed")


## Map it

## Lambert EA: EPSG:3035 ## Lambert CConic: EPSG:3034
tots_sf_prj <- st_as_sf(tots_sf) %>% st_transform(crs = "EPSG:3035")
bbox_prj <- st_bbox( ## manual limits from network data (change as data changes)
  c(xmin = 1271137.2, xmax = 5214780.5, ymax = 5133771.9, ymin = -720268.5),
  crs = "EPSG:3035")

###
map <- ggplot() +
  geom_sf(data = tots_sf_prj, aes(fill = total.ringed), color = "grey20") +
  coord_sf(xlim = 
             c(bbox_prj[1], bbox_prj[3]), ylim = c(bbox_prj[2], bbox_prj[4]), 
           expand = T) +
  scale_fill_continuous(low="#fee0d2", high="#de2d26", 
                        guide="colorbar", na.value="white") +
  # theme_void() +
  theme(
    plot.background = element_rect(fill = "white"),
    # panel.background = element_rect(fill = "slategray1"),
    panel.background = element_rect(fill = "grey85"),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    panel.border = element_rect(colour = "black", fill=NA, size=.5),
    legend.key=element_blank()
  )
map

## SAVE map
ggsave("figures/llimosa_ntl_ringtotalsX.png", plot=map, width=5, height = 6)
