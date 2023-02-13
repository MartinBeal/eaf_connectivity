## Overlay network nodes on cit-sci grid --------------------------------------
# Does number of re-sightings correlate w/ number of checklists?


## project to Europe LAEA (same as cit-sci grid)
node_prj <- st_as_sf(netsf, "nodes") %>% st_transform("EPSG:3035")

node_sv <- vect(node_prj)

x <- terra::extract(arast, node_sv)
node_sv$n_checklist <- x$lyr.1

## filter out first ringing observations (may not be coded...)

node_sv$n_notnew <- node_sv$n_obs - ifelse(is.na(node_sv$n_N), 0, node_sv$n_N)

plot(log(node_sv$n_checklist), log(node_sv$n_notnew))

