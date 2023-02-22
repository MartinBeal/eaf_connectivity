## site 'importance' vs. area/protecteness....

## choose data subset to run --------------------------------------------------

## which network to create
season <- "all"
# season <- "spring"
# season <- "fall"

## Run through each data type -------------------------------------------------
# dtype <- "metal"
# dtype <- "color"
dtype <- "trax"

netsf <- readRDS(
  paste0("data/analysis/networks/", dtype,"_", season, "_iba_hex_10km.rds")
)

## save
allsites <- readRDS(
  "data/analysis/alldatatypes_allsites_poly_all.rds"
)

allnodes <- netsf %>% activate("nodes") %>% sf::st_as_sf()

## add area info to network nodes df ------------------------------------------
allnodes %<>% 
  left_join(st_drop_geometry(allsites[, c("SitRecID", "area", "protect")])) %>% 
  mutate(
    protect = factor( ## doesn't yet include protection level of outsites
      protect, 
    levels = c("little/none", "some", "most", "whole", "unknown", "NA"))
  )

## strength and degree are highly correlated
plot(log(allnodes$degree), log(allnodes$strength))


## site area
plot(log(allnodes$area), (allnodes$strength))
plot(log(allnodes$area), (allnodes$degree))

cor(log(allnodes$area), allnodes$degree)

## site protectedness
ggplot() +
  # geom_boxplot(data=allnodes, aes(x=protect, y=strength, fill=protect)) +
  geom_boxplot(data=allnodes, aes(x=protect, y=degree, fill=protect)) +
  theme_bw() +
  theme(panel.grid = element_blank())

